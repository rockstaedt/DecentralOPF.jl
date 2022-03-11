mutable struct PenaltyTerm
    energy_balance::Vector{Float64}
    upper_flow::Vector{Float64}
    lower_flow::Vector{Float64}
end