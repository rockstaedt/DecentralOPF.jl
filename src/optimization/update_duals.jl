function update_duals!(admm::ADMM)
    __update_lambda!(admm)
    __update_mue!(admm)
    __update_rho!(admm)
end

function __update_lambda!(admm::ADMM)
    lambdas = admm.lambdas[admm.iteration] + admm.gamma * (
        admm.results[admm.iteration].generation
        + admm.results[admm.iteration].discharge
        - admm.results[admm.iteration].charge
        - admm.total_demand
    )
    push!(admm.lambdas, lambdas)
end

function __update_mue!(admm::ADMM)
    injection = admm.results[admm.iteration].injection
    avg_U, _ = get_average_slack_results(admm.iteration)
    mues = (
        admm.mues[admm.iteration]
        + admm.gamma * (admm.ptdf * injection + avg_U .- admm.f_max)
    )
    condition = avg_U .<= 1e-2
    mues .*= condition
    push!(admm.mues, mues)
end

function __update_rho!(admm::ADMM)
    injection = admm.results[admm.iteration].injection
    _, avg_K = get_average_slack_results(admm.iteration)
    rhos = (
        admm.rhos[admm.iteration]
        + admm.gamma * (avg_K- admm.ptdf * injection .- admm.f_max)
    )
    condition = avg_K .<= 1e-2
    rhos .*= condition
    push!(admm.rhos, rhos)
end