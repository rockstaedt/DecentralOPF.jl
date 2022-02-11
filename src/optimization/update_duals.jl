function update_lambda()
    lambdas = admm.lambdas[admm.iteration] + admm.gamma * (
        admm.results[admm.iteration].generation
        + admm.results[admm.iteration].discharge
        - admm.results[admm.iteration].charge
        - admm.total_demand
    )
    push!(admm.lambdas, lambdas)
end

function update_mue()
    injection = admm.results[admm.iteration].injection
    avg_R_ref, _ = get_slack_results(admm.iteration)
    mues = (
        admm.mues[admm.iteration]
        + admm.gamma * (admm.ptdf * injection + avg_R_ref .- admm.L_max)
    )
    condition = avg_R_ref .<= 1e-2
    mues .*= condition
    push!(admm.mues, mues)
end

function update_rho()
    injection = admm.results[admm.iteration].injection
    _, avg_R_cref = get_slack_results(admm.iteration)
    rhos = (
        admm.rhos[admm.iteration]
        + admm.gamma * (avg_R_cref- admm.ptdf * injection .- admm.L_max)
    )
    condition = avg_R_cref .<= 1e-2
    rhos .*= condition
    push!(admm.rhos, rhos)
end

function update_duals()
    update_lambda()
    update_mue()
    update_rho()
end