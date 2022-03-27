function optimize_subproblem(generator::Generator)
    # Create model instance.
    sub = Model(() ->Gurobi.Optimizer(gurobi_env))
    # Omit output in console.
    set_silent(sub)

    # Generation variable with maximum bound.
    @variable(sub, 0 <= P[t=admm.T] <= generator.max_generation)

    # Slack variable for upper flow limit.
    @variable(sub, 0 <= U[line=admm.L, t=admm.T])
    # Slack variable for lower flow limit.
    @variable(sub, 0 <= K[line=admm.L, t=admm.T])

    # Get generator results from previous iteration.
    prev_P = get_unit_results(generator, admm.iteration - 1)

    # Get general node results from previous iteration.
    prev_node_results = get_node_id_to_result(admm.iteration - 1)

    # Define expression for injection term.
    @expression(
        sub,
        injection[n=admm.N, t=admm.T],
        n == admm.node_to_id[generator.node] ? (
            P[t] + (prev_node_results[n].generation[t] - prev_P[t])
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

    # Add all penalty terms to the model.
    add_penalty_terms!(sub)

    # Get current node index.
    node_index = admm.node_to_id[generator.node]

    # Set objective function.
    @objective(
        sub,
        Min,
        sum(P[t] * generator.marginal_costs
            + P[t] * (
                admm.lambdas[admm.iteration][t]
                + sum(
                    admm.ptdf[l, node_index] * (
                        admm.mues[admm.iteration][l, t]
                        - admm.rhos[admm.iteration][l, t]
                    ) for l in admm.L
                )
            )
            + admm.gamma/2 * sub[:penalty_term_eb][t]
            + 10*sub[:penalty_term_upper_flow][t]
            + 10*sub[:penalty_term_lower_flow][t]
            + admm.gamma/2 * sub[:penalty_term_U][t]
            + admm.gamma/2 * sub[:penalty_term_K][t]
            + 1/2 * (P[t] - prev_P[t])^2
            for t in admm.T)
    )

    # Solve optimization problem.
    optimize!(sub)

    # Create penalty term structure.
    penalty_terms = PenaltyTerm(
        value.(sub[:penalty_term_eb]).data,
        value.(sub[:penalty_term_upper_flow]).data,
        value.(sub[:penalty_term_lower_flow]).data
    )

    # Create result structure.
    result = ResultGenerator(
        generator,
        value.(P).data,
        penalty_terms,
        value.(U).data,
        value.(K).data
    )

    return result
end

function optimize_subproblem(storage::Storage)
    # Create model instance.
    sub = Model(() ->Gurobi.Optimizer(gurobi_env))
    # Omit output in console.
    set_silent(sub)

    # Storage variables with maximum bounds.
    @variable(sub, 0 <= D[t=admm.T] <= storage.max_power)
    @variable(sub, 0 <= C[t=admm.T] <= storage.max_power)
    @variable(sub, 0 <= E[t=admm.T] <= storage.max_level)

    # Slack variable for upper flow limit.
    @variable(sub, 0 <= U[line=admm.L, t=admm.T])
    # Slack variable for lower flow limit.
    @variable(sub, 0 <= K[line=admm.L, t=admm.T])

    # Get storage results from previous iteration.
    prev_D, prev_C = get_unit_results(storage, admm.iteration - 1)

    # Get general node results from previous iteration.
    prev_node_results = get_node_id_to_result(admm.iteration - 1)

    # Define expression for injection term.
    @expression(
        sub,
        injection[n=admm.N, t=admm.T],
        n == admm.node_to_id[storage.node] ? (
            prev_node_results[n].generation[t]
            + D[t] + (prev_node_results[n].discharge[t] - prev_D[t])
            - C[t] - (prev_node_results[n].charge[t] - prev_C[t])
            - admm.node_id_to_demand[n][t]
        ) : (
            prev_node_results[n].generation[t]
            + prev_node_results[n].discharge[t]
            - prev_node_results[n].charge[t]
            - admm.node_id_to_demand[n][t]
        )
    )

    # Add all penalty terms to the model.
    add_penalty_terms!(sub)

    # Define constraint for energy level of storage.
    @constraint(
        sub,
        StorageBalance[t=admm.T],
        E[t] == (
            (t == 1 ? 0 : E[t-1]) + C[t] - D[t]
        )
    )

    # Get current node index.
    node_index = admm.node_to_id[storage.node]

    # Set objective function.
    @objective(
        sub,
        Min,
        sum(storage.marginal_costs * (D[t] + C[t])
            + (D[t] - C[t]) * (
                admm.lambdas[admm.iteration][t]
                + sum(
                    admm.ptdf[l, node_index] * (
                        admm.mues[admm.iteration][l, t]
                        - admm.rhos[admm.iteration][l, t]
                    ) for l in admm.L
                )
            )
            + admm.gamma/2 * sub[:penalty_term_eb][t]
            + 10*sub[:penalty_term_upper_flow][t]
            + 10*sub[:penalty_term_lower_flow][t]
            + admm.gamma/2 * sub[:penalty_term_U][t]
            + admm.gamma/2 * sub[:penalty_term_K][t]
            + 1/2 * (D[t] - prev_D[t])^2
            + 1/2 * (C[t] - prev_C[t])^2
            for t in admm.T)
    )

    # Solve optimization problem.
    optimize!(sub)

    # Create penalty term structure.
    penalty_terms = PenaltyTerm(
        value.(sub[:penalty_term_eb]).data,
        value.(sub[:penalty_term_upper_flow]).data,
        value.(sub[:penalty_term_lower_flow]).data
    )

    # Create result structure.
    result = ResultStorage(
        storage,
        value.(D).data,
        value.(C).data,
        value.(E).data,
        penalty_terms,
        value.(U).data,
        value.(K).data
    )

    return result
end