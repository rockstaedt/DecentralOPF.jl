function add_penalty_terms!(sub)
    # Penalty term energy balance
    @expression(
        sub,
        penalty_term_eb[t=admm.T],
        sum(sub[:injection][:, t].data.^2)
    )

    # Penalty term flow1
    @expression(
        sub,
        penalty_term_flow1[t=admm.T],
        sum(
            (
                sum(admm.ptdf[l, n] * sub[:injection][n, t] for n in admm.N)
                + sub[:R_ref][l, t] - admm.L_max[l]
            ).^2 
            for l in admm.L
        )
    )

    # Penalty term flow2
    @expression(
        sub,
        penalty_term_flow2[t=admm.T],
        sum(
            (
                sub[:R_cref][l, t]
                - sum(
                    admm.ptdf[l, n] * sub[:injection].data[n, t]
                    for n in admm.N
                )
                - admm.L_max[l]
            ).^2
            for l in admm.L
        )
    )

    # Penalty terms for slack variables
    avg_R_ref, avg_R_cref = get_slack_results(admm.iteration - 1)

    @expression(
        sub, 
        penalty_term_R_ref[t=admm.T],
        sum((sub[:R_ref].data[:, t] - avg_R_ref[:, t]).^2)
    )

    @expression(
        sub, 
        penalty_term_R_cref[t=admm.T],
        sum((sub[:R_cref].data[:, t] - avg_R_cref[:, t]).^2)
    )
    
end