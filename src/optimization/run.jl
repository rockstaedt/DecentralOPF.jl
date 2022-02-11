function calculate_iteration()
    println("\nIteration: $(admm.iteration)")
    println("###############")
    print_duals(admm.iteration)

    if !isnothing(admm.storages)
        unit_to_result = Dict(zip(
            vcat(generators, storages),
            vcat(
                optimize_subproblem.(generators),
                optimize_subproblem.(storages)
            )
        ))
    else
        unit_to_result = Dict(zip(generators, optimize_subproblem.(generators)))
    end

    result = Result(unit_to_result)

    push!(admm.results, result)

end

function run!(admm)
    while (!admm.convergence.all)
        calculate_iteration()
        
        println("Generation Results: ")
        print_results(true, true, false, admm.iteration)
        
        update_duals()
        
        check_convergence()
    end
end