mutable struct Convergence
    lambda::Bool
    lambda_res::Vector{Vector{Float64}}
    mue::Bool
    mue_res::Vector{Matrix{Float64}}
    rho::Bool
    rho_res::Vector{Matrix{Float64}}
    all::Bool
    function Convergence()
        convergence = new()
        convergence.lambda = false
        convergence.lambda_res = []
        convergence.mue = false
        convergence.mue_res = []
        convergence.rho = false
        convergence.rho_res = []
        return convergence
    end
end