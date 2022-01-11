using JuMP
using Gurobi

# Trick to avoid multiple license printing
gurobi_env = Gurobi.Env()

#-------------------------------------------------------------------------------
# Parameter
#-------------------------------------------------------------------------------

include("example_system_ESR.jl")

D = [180, 370, 70, 570]

T = collect(1:length(D))

#-------------------------------------------------------------------------------
# Main function
#-------------------------------------------------------------------------------

function run_optimization()

end

#-------------------------------------------------------------------------------
# Augmented Lagrangian relaxation with exchange ADMM
#-------------------------------------------------------------------------------

# subproblems 

function optimize_generator_subproblem(generator)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)
    @variable(sub, 0 <= G[t=T] <= gmax[generator])
    @expression(
        sub,
        penalty_term[t=T], 
        (
            G[t] + (
                length(P) * generation_mean[i][t]
                - generation[generator][i-1][t]
            ) 
            + length(S) * discharge_mean[i][t]
            - length(S) * charge_mean[i][t]
            - D[t]
        )^2
    )
    @objective(
        sub,
        Min,
        sum(
            G[t] * mc[g] + lambdas[t]*G[t]
            + gamma / 2 * penalty_term[t] for t in T
        )
    )
    optimize!(sub)
    println("$(s): $(value.(penalty_term).data)")
    return ""
end

function optimize_storage_subproblem(storage)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)
    @variable(sub, 0 <= S_D[t=T] <= gmax[storage])
    @variable(sub, 0 <= S_C[t=T] <= gmax[storage])
    @variable(sub, 0 <= storage_level[t=T] <= capa[storage])
    @expression(
        sub,
        penalty_term[t=T], 
        (
            length(P) * generation_mean[i][t]
            + S_D[t] + (
                length(S) * discharge_mean[i][t]
                - discharge[storage][i-1][t]
            )
            - S_C[t] + (
                length(S) * charge_mean[t]
                - discharge[storage][i-1][t]
            )
            - D[t]
        )^2
    )
    @objective(
        sub,
        Min,
        sum(
            mc[s] * (S_D[t] + S_C[t])
            + lambda[t] * (S_D[t] - S_C[t]) 
            + gamma/2 * penalty_term[t]
            for t in T
        )
    )
    @constraint(
        sub,
        StorageBalance[t=T],
        storage_level[t] == (t == 1 ? 0 : storage_level[t-1]) 
                            + S_C[t] - S_D[t]
    )
    @constraint(
        sub,
        empty_storage,
        storage_level[length(T)] == 0
    )  
    optimize!(sub)
    println("$(s): $(value.(penalty_term).data)")
    return value.(S_D).data, value.(S_C).data, value.(storage_level).data
end

#-------------------------------------------------------------------------------
# Helper Functions
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Optimization
#-------------------------------------------------------------------------------

lambdas = [[3., 3. , 3. , 3.]]
gamma = 0.01

run_optimization()
