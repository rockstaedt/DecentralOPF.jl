using JuMP
using Gurobi

# Trick to avoid multiple license printing
gurobi_env = Gurobi.Env()

include("example_system_ESR.jl")

# Dispatch 

demand = [100, 370, 150, 570]

# Using first element in T as initialization!
T = collect(1:length(demand))

dispatch = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))

@variable(dispatch, 0 <= G[p=P, t=T] <= gmax[p])
@variable(dispatch, 0 <= G_S_d[s=S, t=T] <= gmax[s])
@variable(dispatch, 0 <= G_S_c[s=S, t=T] <= gmax[s])
@variable(dispatch, 0 <= storage_level[s=S, t=T] <= capa[s])

@objective(
    dispatch,
    Min,
    sum(G[p, t] * mc[p] for p in P, t in T)
    + sum(G_S_d[s, t] * mc[s] + G_S_c[s, t] * mc[s] for s in S, t in T)
)

@constraint(
    dispatch,
    EnergyBalance[t=T],
    demand[t] == sum(G[p, t] for p in P) + sum(G_S_d[s, t] for s in S) 
                    - sum(sum(G_S_c[s, t] for s in S)) 
)
@constraint(
    dispatch,
    StorageBalance[s=S, t=T],
    storage_level[s, t] == (t == 1 ? 0 : storage_level[s, t-1])
                            + G_S_c[s, t] - G_S_d[s, t]
)
@constraint(
    dispatch,
    empty_storage[s=S],
    storage_level[s, length(T)] == 0
)

optimize!(dispatch)
objective_value(dispatch)
value.(G)
value.(G_S_d)
value.(G_S_c)
value.(storage_level)

# Augmented Lagrangian relaxation with ADMM

function subproblem_G(g, lambda, G_mean, G_old, G_S_mean)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)
    @variable(sub, 0 <= G[t=T] <= gmax[g])
    @expression(
        sub,
        penalty_term[t=T], 
        (G[t]+(length(P)*G_mean[t]-G_old[t])+length(S)*G_S_mean[t]-demand[t])^2
    )
    @objective(
        sub,
        Min,
        sum(G[t]*mc[g] + lambda[t]*G[t] + gamma/2*penalty_term[t] for t in T)
    )
    optimize!(sub)
    println("$(g): $(value.(penalty_term).data)")
    return value.(G).data
end

function subproblem_G_S(s, lambda, G_S_mean, G_S_old, G_mean)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)
    @variable(sub, -gmax[s] <= G_S[t=T] <= gmax[s])
    @variable(sub, 0 <= storage_level[t=T] <= capa[s])
    @expression(
        sub,
        penalty_term[t=T], 
        (
            length(P) * G_mean[t]
            + G_S[t] + (length(S) * G_S_mean[t] - G_S_old[t])
            -demand[t]
        )^2
    )
    @objective(
        sub,
        Min,
        sum(G_S[t]*mc[s] + lambda[t]*G_S[t]+ gamma/2*penalty_term[t] for t in T)
    )
    @constraint(
        sub,
        StorageBalance[t=T],
        storage_level[t] == (t == 1 ? 0 : storage_level[t-1]) - G_S[t]
    )
    @constraint(
        sub,
        empty_storage,
        storage_level[length(T)-1] == G_S[length(T)]
    )  
    optimize!(sub)
    println("$(s): $(value.(penalty_term).data)")
    return value.(G_S).data, value.(storage_level).data
end

function update_lambda(lambda, G_new, S_new)
    values = []
    for t in T
        push!(
            values,
            lambda[t] + gamma * (
                sum(G_new[p][t] for p in P) + sum(G_S_new[s][t] for s in S)
                - demand[t]
            )
        )
    end
    return values
end

function calc_G_mean(results, amount)
    sum = zeros(Float64, length(T))
    for (g, values) in results
        for n in eachindex(sum)
            sum[n] += values[n]
        end
    end
    return sum ./ amount
end

function check_convergence(lambda)
    epsilon = 10^(-5)
    return sum(abs.(lambda[end]-lambda[length(lambda)-1]).<epsilon) == length(T)
end

begin
    gamma = 0.01
    lambda = [[3., 3., 3., 3.]]
    results = Dict("G" => [], "G_S" => [], "storage_level" => [])
    G_mean = []
    G_S_mean = []
    not_converged = true
    i = 1
    while not_converged
        println("Iteration $(i) with lambda $(lambda[end])")
        global G_new = Dict(
            k => v for (k,v) in zip(P, [
                zeros(Float64, length(T)),
                zeros(Float64, length(T)),
                zeros(Float64, length(T)),
                zeros(Float64, length(T))
            ])
        )
        global G_S_new = Dict(
            k => v for (k,v) in zip(S, [
                zeros(Float64, length(T)),
                zeros(Float64, length(T)),
                zeros(Float64, length(T)),
                zeros(Float64, length(T))
            ])
        )
        global stor_level = Dict(
            k => v for (k,v) in zip(S, [zeros(Float64, length(T))])
        )
        println("Penalty term:")
        for p in P
            if i == 1
                G_new[p] = subproblem_G(
                    p,
                    lambda[end],
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T))
                )
            else
                G_new[p] = subproblem_G(
                    p,
                    lambda[end],
                    G_mean[end],
                    results["G"][end][p],
                    G_S_mean[end]
                )
            end
        end
        for s in S
            if i == 1
                G_S_new[s], stor_level[s] = subproblem_G_S(
                    s,
                    lambda[end],
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T))
                )
            else
                G_S_new[s], stor_level[s] = subproblem_G_S(
                    s,
                    lambda[end],
                    G_S_mean[end],
                    results["G_S"][end][s],
                    G_mean[end]
                )
            end
        end
        println("Generation results:")
        println(G_new)
        push!(results["G"], G_new)
        println("Power storage results:")
        println(G_S_new)
        push!(results["G_S"], G_S_new)
        println("Storage level results:")
        println(stor_level)
        push!(results["storage_level"], stor_level)

        push!(G_mean, calc_G_mean(results["G"][i], length(P)))
        push!(G_S_mean, calc_G_mean(results["G_S"][i], length(S)))

        lambda_new = update_lambda(lambda[i], G_new, G_S_new)
        println("Updated lambda to $(lambda_new)")
        println()
        push!(lambda, lambda_new)
        not_converged = !check_convergence(lambda)
        i += 1
        if i > 10
            break
        end
    end
end
