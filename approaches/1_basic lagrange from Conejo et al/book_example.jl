# Example of Lagrangian Relaxation from book Conejo, p. 196

# min f(x,y) = x^2 + y^2
# st. -x - y <= - 4
#     x >= 0
#     y >= 0
# -> Lagrangian function: L(x,y,lambda) = x^2 + y^2 + lambda*(-x - y + 4)

using JuMP 
using Gurobi

# Trick to avoid multiple license printing
gurobi_env = Gurobi.Env()

function subproblem(variable, lambda)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)
    @variable(sub, variable >= 0);
    @objective(sub, Min, variable^2 - lambda * variable)
    optimize!(sub)
    return value.(variable), objective_value(sub)
end

function updatelambda(lambda, G_new, i)
    a = 1
    b = 0.1
    return (
        lambda 
        + 1/(a+b*i) * (
            (-G_new["x"] - G_new["y"] + 4)
            /(abs(-G_new["x"] - G_new["y"] + 4))
        )
    )
end

function primalproblemvalue(minimizer)
    return minimizer["x"]^2 + minimizer["y"]^2
end

function lagrangianproblemvalue(minimizer, lambda)
    return (
        minimizer["x"]^2 + minimizer["y"]^2 
        + lambda*(-minimizer["x"] - minimizer["y"] + 4)
    )
end


begin
    lambda = [3.]
    not_converged = true
    for i in 1:10
        
        println("Iteration $(i) with lambda $(lambda[end])")
        global G_new = Dict("x" => 0., "y" => 0.)
        global O_new = Dict("x" => 0., "y" => 0.)
        for variable in ["x", "y"]
            G_new[variable], O_new[variable] = subproblem(variable, lambda[end])
        end
        println(G_new)
        println("Primal value: $(primalproblemvalue(G_new))")
        println("Lagrangian value: $(lagrangianproblemvalue(G_new, lambda[end]))")
        
        lambda_new = updatelambda(lambda[end], G_new, i)
        println("Updated lambda to $(lambda_new)")
        push!(lambda, lambda_new)
        println()

    end
end