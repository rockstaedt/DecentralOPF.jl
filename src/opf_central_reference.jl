include("imports.jl")
include("cases/three_node.jl")

# Set up optimization model.
begin
    # Create necessary indexes for all equations.
    N = 1:length(nodes)
    S = 1:length(storages)
    G = 1:length(generators)
    L = 1:length(lines)
    T = 1:length(nodes[1].demand)

    ptdf = calculate_ptdf(nodes, lines)

    # Create model instance.
    m = Model(() ->Gurobi.Optimizer(gurobi_env))

    # Define all optimization variables.

    # Generation power of generators
    @variable(m, 0 <= P[g=G, t=T] <= generators[g].max_generation)
    # Discharge power of storages
    @variable(m, 0 <= D[s=S, t=T] <= storages[s].max_power)
    # Charge power of storages.
    @variable(m, 0 <= C[s=S, t=T] <= storages[s].max_power)
    # Energy level of storages.
    @variable(m, 0 <= E[s=S, t=T] <= storages[s].max_level)
    # Slack variable for upper flow limit.
    @variable(m, 0 <= U[l=L, t=T])
    # Slack variable for lower flow limit.
    @variable(m, 0 <= K[l=L, t=T])
    
    # Define expression for injection term.
    @expression(m, I[n=N, t=T], 
        sum((generators[g].node == nodes[n] ? P[g,t] : 0) for g in G)
        + sum((storages[s].node == nodes[n] ? D[s,t] - C[s,t] : 0) for s in S)
        - nodes[n].demand[t]
    )

    # Set objective function.
    @objective(m, Min, 
        sum(P[g, t] * generators[g].marginal_costs for g in G, t in T) 
        + sum((D[s,t] + C[s,t]) * storages[s].marginal_costs for s in S, t in T) 
    )

    # Set energy balance constraint.
    @constraint(m, EB[t=T], sum(I[:,t]) == 0)
    # Set upper flow constraint.
    @constraint(m, FlowUpper[l=L, t=T], sum(ptdf[l,n] * I[n,t] for n in N) + U[l,t] == lines[l].max_capacity)
    # Set lower flow constraint.
    @constraint(m, FlowLower[l=L, t=T], K[l,t] - sum(ptdf[l,n] * I[n,t] for n in N) == lines[l].max_capacity)
    # Set constraint for energy level of storage.
    @constraint(m, StorageBalance[s=S, t=T], E[s,t] == (t > 1 ? E[s,t-1] : 0) - D[s,t] + C[s,t])
end

# Solve optimization problem.
optimize!(m)

# Get objective value.
println("Objective value: $(objective_value(m))\n")

# Print generation results.
println("Generator results:\n$(value.(P).data)\n")
println("Discharge results:\n$(value.(D).data)\n")
println("Charge results:\n$(value.(C).data)\n")

# Print line flows.
println("Line utilization:\n$(ptdf * value.(I).data)\n")

# Print system price.
lambda = dual.(EB).data
println("System price:\n$(lambda)\n")

# Calculate nodal price.
mu = dual.(FlowUpper).data + dual.(FlowLower).data
nodal_price = zeros(length(N), length(T))
for t in T
    nodal_price[:, t] = lambda[t] .+ sum(mu[l,t] * ptdf[l,:] for l in L)
end
# Print nodal price.
println("Nodal price:\n$(nodal_price)\n")