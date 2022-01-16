using JuMP 
using Gurobi

# Trick to avoid multiple license printing
gurobi_env = Gurobi.Env()

include("example_system.jl")

# Dispatch 

demand = 380

dispatch = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))

@variable(dispatch, 0 <= G[p=P] <= gmax[p])

@objective(dispatch, Min, sum(G[p]^2 * mc[p] + G[p] * mc[p] + mc[p] for p in P));

@constraint(dispatch, EnergyBalance, demand == sum(G[p] for p in P));

optimize!(dispatch)
objective_value(dispatch)
value.(G)

# Lagrangian relaxation

function subproblem(g, lambda)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)
    @variable(sub, 0 <= G <= gmax[g]);
    @objective(sub, Min, G^2 * mc[g] + G * mc[g] - lambda * G)
    optimize!(sub)
    return value.(G), objective_value(sub)
end

function updatelambda(lambda, G_new, i)
    a = 0.01
    b = 0.001
    return (
        lambda 
        + 1/(a+b*i) * (
            (demand - sum(G_new[p] for p in P))
            /(abs(demand - sum(G_new[p] for p in P)))
        )
    )
end

function primalproblemvalue(G)
    return sum(G[p]^2 * mc[p] + G[p] * mc[p] + mc[p] for p in P)
end

function lagrangianproblemvalue(G, lambda)
    return (
        sum(G[p]^2 * mc[p] + G[p] * mc[p] + mc[p] for p in P) 
        + lambda*(demand - sum(G[p] for p in P))
    )
end

# function updatebound(lambda, lower_bound, G_new, O_new)
#     objective_value = (
#         sum(O_new[p] for p in P) + lambda * (demand - sum(G_new[p] for p in P))
#     )
#     if objective_value > lower_bound
#         return objective_value
#     else
#         return lower_bound
#     end
# end

begin
    lambda = [6000.]
    for i in 1:400
        
        println("Iteration $(i) with lambda $(lambda[end])")
        global G_new = Dict(k => v for (k,v) in zip(P, [0. 0. 0. 0.]))
        global O_new = Dict(k => v for (k,v) in zip(P, [0. 0. 0. 0.]))
        for p in P
            G_new[p], O_new[p] = subproblem(p, lambda[end])
        end
        println(G_new)
        println("Primal value: $(primalproblemvalue(G_new))")
        println("Lagrangian value: $(lagrangianproblemvalue(G_new, lambda[end]))")

        lambda_new = updatelambda(lambda[end], G_new, i)
        println("Lambda diff: $(abs(lambda_new-lambda[end]))")
        println("Updated lambda to $(lambda_new)")
        println()

        push!(lambda, lambda_new)
        if sum(G_new[p] for p in P) > demand
            break
        end
    end
end
