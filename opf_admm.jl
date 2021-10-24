using JuMP
using Gurobi

# Trick to avoid multiple license printing
gurobi_env = Gurobi.Env()

include("example_system.jl")

# Dispatch 

demand = 380

dispatch = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))

@variable(dispatch, 0 <= G[p=P] <= gmax[p])

@objective(dispatch, Min, sum(G[p] * mc[p] for p in P));

@constraint(dispatch, EnergyBalance, demand == sum(G[p] for p in P));

optimize!(dispatch)
objective_value(dispatch)
value.(G)

# Augmented Lagrangian relaxation with ADMM

function subproblem(g, lambda, G_mean, G_old)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)
    @variable(sub, 0 <= G <= gmax[g]);
    @expression(sub, penalty_term, 
        (G + (G_mean - G_old) - demand)^2
    );
    @objective(sub, Min, G * mc[g] + lambda * G + gamma/2 * penalty_term)
    optimize!(sub)
    return value.(G)
end

function update_lambda(lambda, G_new)
    return lambda + gamma * (sum(G_new[p] for p in P) - demand)
end

function calc_G_mean(i)
    sum = 0
    for (key, value) in results[i]
        sum = sum + value
    end
    return sum / length(P)
end

function check_convergence(lambda)
    epsilon = 10^(-5)
    return abs(lambda[end] - lambda[length(lambda) - 1]) < epsilon
end

begin
    gamma = 0.3
    lambda = [3.]
    results = []
    G_mean = []
    not_converged = true
    i = 1
    while not_converged
        println("Iteration $(i) with lambda $(lambda[end])")
        G_new = Dict(k => v for (k,v) in zip(P, [0., 0., 0., 0.]))
        for p in P
            if i == 1
                G_new[p] = subproblem(p, lambda[end], 0, 0)
            else
                G_new[p] = subproblem(
                    p,
                    lambda[end],
                    G_mean[end],
                    results[end][p]
                )

            end
        end
        println(G_new)
        push!(results, G_new)
        push!(G_mean, calc_G_mean(i))

        lambda_new = update_lambda(lambda[end], G_new)
        println("Lambda diff: $(abs(lambda_new-lambda[end]))")
        println("Updated lambda to $(lambda_new)")
        println()
        push!(lambda, lambda_new)
        not_converged = !check_convergence(lambda)
        i += 1
    end
end
