function add_penalty_terms!(sub)
    # Penalty term energy balance.
    @expression(
        sub,
        penalty_term_eb[t=admm.T],
        sum(sub[:injection][:, t].data)^2
    )

    # Penalty term upper flow.
    @expression(
        sub,
        penalty_term_upper_flow[t=admm.T],
        sum(
            (
                sum(admm.ptdf[l, n] * sub[:injection][n, t] for n in admm.N)
                + sub[:U][l, t] - admm.f_max[l]
            ).^2 
            for l in admm.L
        )
    )

    # Penalty term lower flow.
    @expression(
        sub,
        penalty_term_lower_flow[t=admm.T],
        sum(
            (
                sub[:K][l, t]
                - sum(
                    admm.ptdf[l, n] * sub[:injection].data[n, t]
                    for n in admm.N
                )
                - admm.f_max[l]
            ).^2
            for l in admm.L
        )
    )

    # Penalty terms for slack variables.
    avg_U, avg_K = get_average_slack_results(admm.iteration - 1)

    @expression(
        sub, 
        penalty_term_U[t=admm.T],
        sum((sub[:U].data[:, t] - avg_U[:, t]).^2)
    )

    @expression(
        sub, 
        penalty_term_K[t=admm.T],
        sum((sub[:K].data[:, t] - avg_K[:, t]).^2)
    )
    
end