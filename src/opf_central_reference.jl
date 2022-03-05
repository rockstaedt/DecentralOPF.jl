include("imports.jl")
include("cases/three_node.jl")

begin # Like ADMM
    N = 1:length(nodes)
    S = 1:length(storages)
    P = 1:length(generators)
    L = 1:length(lines)
    T = 1:length(nodes[1].demand)

    ptdf = calculate_ptdf(nodes, lines)

    m = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    @variable(m, 0 <= G[t=T, p=P] <= generators[p].max_generation)
    @variable(m, 0 <= G_S_d[t=T, s=S] <= storages[s].max_power)
    @variable(m, 0 <= G_S_c[t=T, s=S] <= storages[s].max_power)
    @variable(m, 0 <= L_S[t=T, s=S] <= storages[s].max_level)
    @variable(m, 0 <= R_ref[t=T, l=L])
    @variable(m, 0 <= R_cref[t=T, l=L])
    
    @expression(m, I[t=T, n=N], 
        sum((generators[p].node == nodes[n] ? G[t,p] : 0) for p in P)
        + sum((storages[s].node == nodes[n] ? G_S_d[t,s] - G_S_c[t,s] : 0) for s in S)
        - nodes[n].demand[t]
    )
    @objective(m, Min, 
        sum(G[t,p] * generators[p].marginal_costs for p in P, t in T) 
        + sum((G_S_d[t,s] + G_S_c[t,s]) * storages[s].marginal_costs for s in S, t in T) 
    )
    @constraint(m, EB[t=T], sum(I[t, :]) == 0)
    @constraint(m, FlowUpper[t=T, l=L], ptdf[l,:]' * I[t, :].data + R_ref[t, l] == lines[l].max_capacity)
    @constraint(m, FlowLower[t=T, l=L], R_cref[t, l] - ptdf[l,:]' * I[t, :].data == lines[l].max_capacity)
    @constraint(m, StorageBalance[t=T, s=S], L_S[t,s] == (t > 1 ? L_S[t-1, s] : 0) - G_S_d[t,s] + G_S_c[t,s])
end

optimize!(m)
objective_value(m)
print(value.(G))
value.(G_S_d)
value.(G_S_c)
print(ptdf * value.(I).data')
lambda = dual.(EB).data
mu = dual.(FlowUpper).data + dual.(FlowLower).data
nodal_price = zeros(length(N), length(T))
for t in T
    nodal_price[:, t] = lambda[t] .+ sum(mu[t, l] * ptdf[l, :] for l in L)
end
println(nodal_price)