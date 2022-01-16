using JuMP
using Gurobi
using PlotlyJS

include("structures.jl")


# Trick to avoid multiple license printing
gurobi_env = Gurobi.Env()
begin 
    node1 = Node("N1", [10, 250], false)
    node2 = Node("N2", [50, 30], false)
    node3 = Node("N3", [50, 150], true)

    nodes = [node1, node2, node3]

    line1 = Line("L1", node2, node1, 30, 1)
    line2 = Line("L2", node3, node1, 100, 1)
    line3 = Line("L3", node2, node3, 100, 2)

    lines = [line1, line2, line3]

    pv = Generator("pv", 3, 80, "yellow", node1)
    wind = Generator("wind", 4, 120, "lightblue", node2)
    coal = Generator("coal", 30, 300, "brown", node3)
    gas = Generator("gas", 50, 120, "grey", node1)

    generators = [pv, wind, coal, gas]

    battery = Storage("battery", 1, 10, 20, "purple", node1)
    ps = Storage("ps", 2, 15, 30, "blue", node1)

    storages = [battery]
end

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
    @constraint(m, FlowUpper[t=T, l=L], ptdf[l,:]' * I[t, :].data <= lines[l].max_capacity)
    @constraint(m, FlowLower[t=T, l=L], ptdf[l,:]' * I[t, :].data >= -lines[l].max_capacity)
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

begin # Like normal OPF with lambda[t,n] and variable Injection
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
    @variable(m, I[T, N])
    
    @objective(m, Min, 
        sum(G[t,p] * generators[p].marginal_costs for p in P, t in T) 
        + sum((G_S_d[t,s] + G_S_c[t,s]) * storages[s].marginal_costs for s in S, t in T) 
    )

    @constraint(m, EB[t=T, n=N], 
        sum((generators[p].node == nodes[n] ? G[t,p] : 0) for p in P)
        + sum((storages[s].node == nodes[n] ? G_S_d[t,s] - G_S_c[t,s] : 0) for s in S)
        - nodes[n].demand[t] == I[t,n]
    )
    @constraint(m, [t=T], sum(I[t, :]) == 0)
    @constraint(m, FlowUpper[t=T, l=L], ptdf[l,:]' * I[t, :].data <= lines[l].max_capacity)
    @constraint(m, FlowLower[t=T, l=L], ptdf[l,:]' * I[t, :].data >= -lines[l].max_capacity)
    @constraint(m, StorageBalance[t=T, s=S], L_S[t,s] == (t > 1 ? L_S[t-1, s] : 0) - G_S_d[t,s] + G_S_c[t,s])
end


optimize!(m)
objective_value(m)
value.(G)
value.(G_S_d)
value.(G_S_c)
ptdf * value.(I).data'

lambda = dual.(EB).data'
lambda
nodal_price
