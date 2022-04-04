function run!(admm::ADMM)
    while (!admm.convergence.all)
        calculate_iteration!(admm)
    end
end

function calculate_iteration!(admm::ADMM)
    optimize_all_subproblems!(admm)

    println("Generation Results: ")
    print_results(true, true, false, admm.iteration)
    
    update_duals!(admm)
    
    check_convergence!(admm)
end