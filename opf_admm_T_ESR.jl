using JuMP
using Gurobi

# Trick to avoid multiple license printing
gurobi_env = Gurobi.Env()

include("example_system_ESR.jl")

# Dispatch 

demand = [50, 80, 250, 600]

# Using first element in T as initialization!
T = collect(1:length(demand))

dispatch = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))

@variable(dispatch, 0 <= G[p=P, t=T] <= gmax[p])
@variable(dispatch, -gmax[s] <= G_S[s=S, t=T] <= gmax[s])
@variable(dispatch, 0 <= storage_level[s=S, t=T] <= capa[s])

@objective(
    dispatch,
    Min,
    sum(G[p, t] * mc[p] for p in P, t in T)
    + sum(G_S[s, t] * mc[s] for s in S, t in T)
)

@constraint(
    dispatch,
    EnergyBalance[t=T],
    demand[t] == sum(G[p, t] for p in P) + sum(G_S[s, t] for s in S)
)
@constraint(
    dispatch,
    StorageBalance[s=S, t=T],
    storage_level[s, t] == (t == 1 ? 0 : storage_level[s, t-1]) - G_S[s, t]
)
@constraint(
    dispatch,
    empty_storage[s=S],
    storage_level[s,
    length(T)-1] == G_S[s, length(T)]
)

optimize!(dispatch)
objective_value(dispatch)
value.(G)
value.(G_S)
value.(storage_level)

# Augmented Lagrangian relaxation with ADMM

# function subproblem(g, lambda, G_mean, G_old)
#     sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
#     set_silent(sub)
#     @variable(sub, 0 <= G[t=T] <= gmax[g])
#     @objective(
#         sub,
#         Min,
#         sum(
#             G[t] * mc[g] + lambda[t] * G[t]
#             + gamma/2 * ((G[t] + (G_mean[t] - G_old[t]) - demand[t])^2)
#             for t in T
#         )
#     )
#     optimize!(sub)
#     return value.(G).data
# end

# function update_lambda(lambda, G_new)
#     values = []
#     for t in T
#         push!(
#             values,
#             lambda[t] + gamma * (sum(G_new[p][t] for p in P) - demand[t])
#         )
#     end
#     return values
# end

# function calc_G_mean(i)
#     sum = zeros(Float64, length(T))
#     for (g, values) in results[i]
#         for n in eachindex(sum)
#             sum[n] += values[n]
#         end
#     end
#     return sum ./ length(P)
# end

# function check_convergence(lambda)
#     epsilon = 10^(-5)
#     return sum(abs.(lambda[end]-lambda[length(lambda)-1]).<epsilon) == length(T)
# end

# begin
#     gamma = 0.3
#     lambda = [[3., 3., 3., 3.]]
#     results = []
#     G_mean = []
#     not_converged = true
#     i = 1
#     while not_converged
#         println("Iteration $(i) with lambda $(lambda[end])")
#         global G_new = Dict(
#             k => v for (k,v) in zip(P, [
#                 zeros(Float64, length(T)),
#                 zeros(Float64, length(T)),
#                 zeros(Float64, length(T)),
#                 zeros(Float64, length(T))
#             ])
#         )
#         for p in P
#             if i == 1
#                 G_new[p] = subproblem(p, lambda[end], [0,0,0,0], [0,0,0,0])
#             else
#                 G_new[p] = subproblem(
#                     p,
#                     lambda[end],
#                     G_mean[end],
#                     results[end][p]
#                 )
#             end
#         end
#         println(G_new)
#         push!(results, G_new)
#         push!(G_mean, calc_G_mean(i))

#         lambda_new = update_lambda(lambda[i], G_new)
#         println("Updated lambda to $(lambda_new)")
#         println()
#         push!(lambda, lambda_new)
#         not_converged = !check_convergence(lambda)
#         i += 1
#     end
# end
