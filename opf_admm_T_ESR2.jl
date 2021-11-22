using JuMP
using Gurobi
using Plots

# Trick to avoid multiple license printing
gurobi_env = Gurobi.Env()

include("example_system_ESR.jl")

# Dispatch 

D = [160, 220, 140, 240]

# Using first element in T as initialization!
T = collect(1:length(D))

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
    D[t] == sum(G[p, t] for p in P) + sum(G_S_d[s, t] for s in S) 
                    - sum(G_S_c[s, t] for s in S)
)
@constraint(
    dispatch,
    StorageBalance[s=S, t=T],
    storage_level[s, t] == (t == 1 ? 0 : storage_level[s, t-1])
                            + G_S_c[s, t] - G_S_d[s, t]
)
# @constraint(
#     dispatch,
#     empty_storage[s=S],
#     storage_level[s, length(T)] == 0
# )

optimize!(dispatch)
objective_value(dispatch)
value.(G)
value.(G_S_d)
value.(G_S_c)
value.(storage_level)
dual.(EnergyBalance)

# Augmented Lagrangian relaxation with ADMM

function subproblem_G(g, lambda, G_old, G_S_d_old, G_S_c_old)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)
    @variable(sub, 0 <= G[t=T] <= gmax[g])
    @expression(
        sub,
        penalty_term[t=T], 
        (
            G[t] + sum(G_old[p][t] for p in setdiff(P, [g])) 
            + sum(G_S_d_old[s][t] for s in S)
            - sum(G_S_c_old[s][t] for s in S)
            - D[t]
        )^2
    )
    @objective(
        sub,
        Min,
        sum(G[t]*mc[g] + lambda[t]*G[t] + gamma/2*penalty_term[t] for t in T)
    )
    optimize!(sub)
    println("$(g): $(value.(penalty_term).data)")
    return value.(G).data, value.(penalty_term).data
end

function subproblem_G_S(s,lambda, G_S_d_old, G_S_c_old, G_old)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)
    @variable(sub, 0 <= G_S_d[t=T] <= gmax[s])
    @variable(sub, 0 <= G_S_c[t=T] <= gmax[s])
    @variable(sub, 0 <= storage_level[t=T] <= capa[s])
    @expression(
        sub,
        penalty_term[t=T], 
        (
            sum(G_old[p][t] for p in P)
            + G_S_d[t] + sum(G_S_d_old[stor][t] for stor in setdiff(S, [s]))
            - G_S_c[t] - sum(G_S_c_old[stor][t] for stor in setdiff(S, [s]))
            - D[t]
        )^2
    )
    @objective(
        sub,
        Min,
        sum(
            mc[s] * (G_S_d[t] + G_S_c[t])
            + lambda[t] * (G_S_d[t] - G_S_c[t]) 
            + gamma/2 * penalty_term[t]
            for t in T
        )
    )
    @constraint(
        sub,
        StorageBalance[t=T],
        storage_level[t] == (t == 1 ? 0 : storage_level[t-1]) 
                            + G_S_c[t] - G_S_d[t]
    )
    # @constraint(
    #     sub,
    #     empty_storage,
    #     storage_level[length(T)] == 0
    # )  
    optimize!(sub)
    println("$(s): $(value.(penalty_term).data)")
    return value.(G_S_d).data, value.(G_S_c).data, value.(storage_level).data, value.(penalty_term).data
end

function update_lambda(lambda, G_new, G_S_d_new, G_S_c_new)
    values = []
    for t in T
        push!(
            values,
            lambda[t] + gamma * (
                sum(G_new[p][t] for p in P) + sum(G_S_d_new[s][t] for s in S)
                - sum(G_S_c_new[s][t] for s in S)
                - D[t]
            )
        )
    end
    return values
end

function check_convergence(lambda)
    epsilon = 10^(-5)
    return sum(abs.(lambda[end]-lambda[length(lambda)-1]).<epsilon) == length(T)
end

begin
    gamma = 0.08
    lambda = [[3., 3., 3., 3.]]
    penalty_term = []
    results = Dict(
        "G" => [
            Dict(
                k => v for (k,v) in zip(P, [
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T))
                ])
            )
        ],
        "G_S_d" => [
            Dict(
                k => v for (k,v) in zip(S, [
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T))
                ])
            )
        ],
        "G_S_c" => [
            Dict(
                k => v for (k,v) in zip(S, [
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T)),
                    zeros(Float64, length(T))
                ])
            )
        ],
        "storage_level" => []
    )
    not_converged = true
    i = 1
    while true
        println("Iteration $(i) with lambda $(lambda[end])")
        penalties_unit = []
        global G_new = Dict(
            k => v for (k,v) in zip(P, [
                zeros(Float64, length(T)),
                zeros(Float64, length(T)),
                zeros(Float64, length(T)),
                zeros(Float64, length(T))
            ])
        )
        global G_S_d_new = Dict(
            k => v for (k,v) in zip(S, [
                zeros(Float64, length(T)),
                zeros(Float64, length(T)),
                zeros(Float64, length(T)),
                zeros(Float64, length(T))
            ])
        )
        global G_S_c_new = Dict(
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
            G_new[p], penalties = subproblem_G(
                p,
                lambda[end],
                results["G"][end],
                results["G_S_d"][end],
                results["G_S_c"][end]
            )
            push!(penalties_unit, penalties)
        end
        println("Generation results:")
        println(G_new)
        push!(results["G"], G_new)
        for s in S
            G_S_d_new[s], G_S_c_new[s], stor_level[s], penalties = subproblem_G_S(
                s,
                lambda[end],
                results["G_S_d"][end],
                results["G_S_c"][end],
                results["G"][end]
            )
            push!(penalties_unit, penalties)
        end
        penalty_sum = [0., 0., 0., 0.]
        for penalties in penalties_unit
            for t in T
                penalty_sum[t] += penalties[t]
            end
        end
        push!(penalty_term, penalty_sum)
        println("Storage discharge results:")
        println(G_S_d_new)
        push!(results["G_S_d"], G_S_d_new)
        println("Storage charge results:")
        println(G_S_c_new)
        push!(results["G_S_c"], G_S_c_new)
        println("Storage level results:")
        println(stor_level)
        push!(results["storage_level"], stor_level)

        lambda_new = update_lambda(lambda[i], G_new, G_S_d_new, G_S_c_new)
        # if lambda_new[4] < -30.
        #     gamma = 0.025
        # end
        println("Updated lambda to $(lambda_new)")
        println()
        push!(lambda, lambda_new)
        if !not_converged
            break
        end
        not_converged = !check_convergence(lambda)
        i += 1
        if i > 800
            break
        end
    end
end

result_type = "G_S_c"
unit = "esr"
x = []
y = []
for t in T
    push!(x, 1:i-1)
    push!(y, [])
end
for iteration in 1:i-1
    for t in T
        push!(y[t], results[result_type][iteration][unit][t])
    end
end
plot(x, y, labels=["t1" "t2" "t3" "t4"], titles=result_type * " - " * unit)

result_type = "G_S_d"
unit = "esr"
x = []
y = []
for t in T
    push!(x, 1:i-1)
    push!(y, [])
end
for iteration in 1:i-1
    for t in T
        push!(y[t], results[result_type][iteration][unit][t])
    end
end
plot(x, y, labels=["t1" "t2" "t3" "t4"], titles=result_type * " - " * unit)

x = []
y = []
for t in T
    push!(x, 1:i-1)
    push!(y, [])
end
for iteration in 1:i-1
    for t in T
        push!(y[t], lambda[iteration][t])
    end
end
plot(x, y, labels=["t1" "t2" "t3" "t4"], titles="lambda")

x = []
y = []
for t in T
    push!(x, 1:i-1)
    push!(y, [])
end
for iteration in 1:i-1
    for t in T
        push!(y[t], penalty_term[iteration][t])
    end
end
plot(x, y, labels=["t1" "t2" "t3" "t4"], titles="penalty_term")