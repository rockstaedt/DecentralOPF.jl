function sum_up(summand1::PenaltyTerm, summand2::PenaltyTerm)
    summand1.energy_balance += summand2.energy_balance
    summand1.upper_flow += summand2.upper_flow
    summand1.lower_flow += summand2.lower_flow
    return summand1
end

function get_empty_penalty_term()
    return PenaltyTerm(
        zeros(Float64, length(admm.T)),
        zeros(Float64, length(admm.T)),
        zeros(Float64, length(admm.T))
    )
end