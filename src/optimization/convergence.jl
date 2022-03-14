function check_convergence!(admm::ADMM)
    eps = 10^(-3)
    if admm.iteration != 1
        lambda_res = abs.(
            (
                admm.results[admm.iteration].generation
                + admm.results[admm.iteration].discharge
                - admm.results[admm.iteration].charge
                - admm.total_demand
            )
        )
        push!(admm.convergence.lambda_res, lambda_res)

        mue_res = abs.(
            (
                admm.ptdf * admm.results[admm.iteration].injection
                + admm.results[admm.iteration].avg_U
                .- admm.f_max
            )
        )
        push!(admm.convergence.mue_res, mue_res)

        rho_res = abs.(
            (
                admm.results[admm.iteration].avg_K
                - admm.ptdf * admm.results[admm.iteration].injection
                .- admm.f_max
            )
        )
        push!(admm.convergence.rho_res, rho_res)

        admm.convergence.lambda = all(lambda_res .< eps)
        admm.convergence.mue = all(mue_res .< eps)
        admm.convergence.rho = all(rho_res .< eps)
        
        admm.convergence.all = all([
            admm.convergence.lambda
            admm.convergence.mue
            admm.convergence.rho
        ])
    end
    
    if admm.convergence.all
        println("Converged")
    else
        println("Not converged")
        admm.iteration += 1
    end
end