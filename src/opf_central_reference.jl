include("imports.jl")
include("cases/three_node.jl")

# Set up optimization model.
begin
    # Create necessary indexes for all equations.
    N = 1:length(nodes)
    S = 1:length(storages)
    P = 1:length(generators)
    L = 1:length(lines)
    T = 1:length(nodes[1].demand)

    ptdf = calculate_ptdf(nodes, lines)

    # Create model instance.
    m = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))

    # Define all optimization variables.

    # Generation power of generators
    @variable(m, 0 <= G[p=P, t=T] <= generators[p].max_generation)
    # Discharge power of storages
    @variable(m, 0 <= G_S_d[s=S, t=T] <= storages[s].max_power)
    # Charge power of storages.
    @variable(m, 0 <= G_S_c[s=S, t=T] <= storages[s].max_power)
    # Energy level of storages.
    @variable(m, 0 <= L_S[s=S, t=T] <= storages[s].max_level)
    # Slack variable for upper flow limit.
    @variable(m, 0 <= R_ref[l=L, t=T])
    # Slack variable for lower flow limit.
    @variable(m, 0 <= R_cref[l=L, t=T])
    
    # Define expression for injection term.
    @expression(m, I[n=N, t=T], 
        sum((generators[p].node == nodes[n] ? G[p,t] : 0) for p in P)
        + sum((storages[s].node == nodes[n] ? G_S_d[s,t] - G_S_c[s,t] : 0) for s in S)
        - nodes[n].demand[t]
    )

    # Set objective function.
    @objective(m, Min, 
        sum(G[p, t] * generators[p].marginal_costs for p in P, t in T) 
        + sum((G_S_d[s,t] + G_S_c[s,t]) * storages[s].marginal_costs for s in S, t in T) 
    )

    # Set energy balance constraint.
    @constraint(m, EB[t=T], sum(I[:,t]) == 0)
    # Set upper flow constraint.
    @constraint(m, FlowUpper[l=L, t=T], sum(ptdf[l,n] * I[n,t] for n in N) + R_ref[l,t] == lines[l].max_capacity)
    # Set lower flow constraint.
    @constraint(m, FlowLower[l=L, t=T], R_cref[l,t] - sum(ptdf[l,n] * I[n,t] for n in N) == lines[l].max_capacity)
    # Set constraint for energy level of storage.
    @constraint(m, StorageBalance[s=S, t=T], L_S[s,t] == (t > 1 ? L_S[s,t-1] : 0) - G_S_d[s,t] + G_S_c[s,t])
end

# Solve optimization problem.
optimize!(m)

# Get objective value.
println("Objective value: $(objective_value(m))\n")

# Print generation results.
println("Generator results:\n$(value.(G)')\n")
println("Discharge results:\n$(value.(G_S_d)')\n")
println("Charge results:\n$(value.(G_S_c)')\n")

# Print line flows.
println("Line utilization:\n$(ptdf * value.(I).data')\n")

# Print system price.
lambda = dual.(EB).data
println("System price:\n$(lambda)\n")

# Calculate nodal price.
mu = dual.(FlowUpper).data + dual.(FlowLower).data
nodal_price = zeros(length(N), length(T))
for t in T
    nodal_price[:, t] = lambda[t] .+ sum(mu[l,t] * ptdf[l,:] for l in L)
end
println("Nodal price:\n$(nodal_price)\n")