function optimize_subproblem(generator::Generator)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)

    # Maximum capacity of generator
    @variable(sub, 0 <= G[t=admm.T] <= generator.max_generation)

    # Slack variables
    @variable(sub, 0 <= R_ref[line=admm.L, t=admm.T])
    @variable(sub, 0 <= R_cref[line=admm.L, t=admm.T])


    prev_G = get_unit_results(generator, admm.iteration - 1)

    prev_node_results = get_node_id_to_result(admm.iteration - 1)

    # Injection term
    @expression(
        sub,
        injection[n=admm.N, t=admm.T],
        n == admm.node_to_id[generator.node] ? (
            G[t] + (prev_node_results[n].generation[t] - prev_G[t])
            + prev_node_results[n].discharge[t]
            - prev_node_results[n].charge[t]
            - admm.node_id_to_demand[n][t]
        ) : (
            prev_node_results[n].generation[t]
            + prev_node_results[n].discharge[t]
            - prev_node_results[n].charge[t]
            - admm.node_id_to_demand[n][t]
        )
    )

    add_penalty_terms!(sub)

    node_index = admm.node_to_id[generator.node]

    @objective(
        sub,
        Min,
        sum(G[t] * generator.marginal_costs
            + G[t] * (
                admm.lambdas[admm.iteration][t]
                + sum(
                    admm.ptdf[l, node_index] * (
                        admm.mues[admm.iteration][l, t]
                        - admm.rhos[admm.iteration][l, t]
                    ) for l in admm.L
                )
            )
            + admm.gamma/2 * sub[:penalty_term_eb][t]
            + 10*sub[:penalty_term_flow1][t]
            + 10*sub[:penalty_term_flow2][t]
            + admm.gamma/2 * sub[:penalty_term_R_ref][t]
            + admm.gamma/2 * sub[:penalty_term_R_cref][t]
            + 1/2 * (G[t] - prev_G[t])^2
            for t in admm.T)
    )

    optimize!(sub)

    penalty_terms = PenaltyTerm(
        value.(sub[:penalty_term_eb]).data,
        value.(sub[:penalty_term_flow1]).data,
        value.(sub[:penalty_term_flow2]).data
    )

    result = ResultGenerator(
        generator,
        value.(G).data,
        penalty_terms,
        value.(R_ref).data,
        value.(R_cref).data
    )

    return result
end

function optimize_subproblem(storage::Storage)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)

    # Storage variables with maximum bounds
    @variable(sub, 0 <= G_S_d[t=admm.T] <= storage.max_power)
    @variable(sub, 0 <= G_S_c[t=admm.T] <= storage.max_power)
    @variable(sub, 0 <= storage_level[t=admm.T] <= storage.max_level)

    # Slack variables
    @variable(sub, 0 <= R_ref[line=admm.L, t=admm.T])
    @variable(sub, 0 <= R_cref[line=admm.L, t=admm.T])


    previous_S_d, previous_S_c = get_unit_results(storage, admm.iteration - 1)

    prev_node_results = get_node_id_to_result(admm.iteration - 1)

    # Injection term
    @expression(
        sub,
        injection[n=admm.N, t=admm.T],
        n == admm.node_to_id[storage.node] ? (
            prev_node_results[n].generation[t]
            + G_S_d[t] + (prev_node_results[n].discharge[t] - previous_S_d[t])
            - G_S_c[t] - (prev_node_results[n].charge[t] - previous_S_c[t])
            - admm.node_id_to_demand[n][t]
        ) : (
            prev_node_results[n].generation[t]
            + prev_node_results[n].discharge[t]
            - prev_node_results[n].charge[t]
            - admm.node_id_to_demand[n][t]
        )
    )

    add_penalty_terms!(sub)

    # Storage level constraint
    @constraint(
        sub,
        StorageBalance[t=admm.T],
        storage_level[t] == (
            (t == 1 ? 0 : storage_level[t-1]) + G_S_c[t] - G_S_d[t]
        )
    )

    node_index = admm.node_to_id[storage.node]

    @objective(
        sub,
        Min,
        sum(storage.marginal_costs * (G_S_d[t] + G_S_c[t])
            + (G_S_d[t] - G_S_c[t]) * (
                admm.lambdas[admm.iteration][t]
                + sum(
                    admm.ptdf[l, node_index] * (
                        admm.mues[admm.iteration][l, t]
                        - admm.rhos[admm.iteration][l, t]
                    ) for l in admm.L
                )
            )
            + admm.gamma/2 * sub[:penalty_term_eb][t]
            + 10*sub[:penalty_term_flow1][t]
            + 10*sub[:penalty_term_flow2][t]
            + admm.gamma/2 * sub[:penalty_term_R_ref][t]
            + admm.gamma/2 * sub[:penalty_term_R_cref][t]
            + 1/2 * (G_S_d[t] - previous_S_d[t])^2
            + 1/2 * (G_S_c[t] - previous_S_c[t])^2
            for t in admm.T)
    )

    optimize!(sub)

    penalty_terms = PenaltyTerm(
        value.(sub[:penalty_term_eb]).data,
        value.(sub[:penalty_term_flow1]).data,
        value.(sub[:penalty_term_flow2]).data
    )

    result = ResultStorage(
        storage,
        value.(G_S_d).data,
        value.(G_S_c).data,
        value.(storage_level).data,
        penalty_terms,
        value.(R_ref).data,
        value.(R_cref).data
    )

    return result
end