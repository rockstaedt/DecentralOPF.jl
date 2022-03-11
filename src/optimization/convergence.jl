function check_convergence()
    eps = 10^(-3)
    if admm.iteration != 1
        r_delta = abs.(
            (
                admm.results[admm.iteration].generation
                + admm.results[admm.iteration].discharge
                - admm.results[admm.iteration].charge
                - admm.total_demand
            )
        )
        push!(admm.convergence.lambda_res, r_delta)
        s_delta = abs.(
            (
                admm.ptdf * admm.results[admm.iteration].injection
                + admm.results[admm.iteration].avg_U
                .- admm.f_max
            )
        )
        push!(admm.convergence.mue_res, s_delta)
        t_delta = abs.(
            (
                admm.results[admm.iteration].avg_K
                - admm.ptdf * admm.results[admm.iteration].injection
                .- admm.f_max
            )
        )
        push!(admm.convergence.rho_res, t_delta)
        admm.convergence.lambda = all(r_delta .< eps)
        admm.convergence.mue = all(s_delta .< eps)
        admm.convergence.rho = all(t_delta .< eps)
        
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