using DataFrames

mutable struct Generator
    name::String
    marginal_costs::Int
    max_generation::Int
    plot_color::String
end

mutable struct Storage
    name::String
    marginal_costs::Int
    max_power::Int
    max_level::Int
    plot_color::String
end

struct ResultStorage
    storage::Storage
    discharge::Vector{Float64}
    charge::Vector{Float64}
    level::Vector{Float64}
    penalty_term::Vector{Float64}
end

struct ResultGenerator
    generator::Generator
    generation::Vector{Float64}
    penalty_term::Vector{Float64}
end

mutable struct Result
    unit_to_result::Dict
    sum_G_t::Vector{Float64}
    sum_S_d_t::Vector{Float64}
    sum_S_c_t::Vector{Float64}
    total_costs::Float64
    sum_penalty_t::Vector{Float64}
    
    function Result(unit_to_result::Dict)
        result = new()

        result.unit_to_result = unit_to_result
        result.sum_G_t = zeros(Float64, length(admm.demand))
        result.sum_S_d_t = zeros(Float64, length(admm.demand))
        result.sum_S_c_t = zeros(Float64, length(admm.demand))
        result.total_costs = 0
        result.sum_penalty_t = zeros(Float64, length(admm.demand))

        for (unit, result_unit) in unit_to_result
            result.sum_penalty_t += result_unit.penalty_term
            if typeof(unit) == Storage
                result.sum_S_d_t += result_unit.discharge
                result.sum_S_c_t += result_unit.charge
                result.total_costs += result_unit.storage.marginal_costs * (
                    sum(result_unit.discharge + result_unit.charge)
                )
            elseif typeof(unit) == Generator
                result.sum_G_t += result_unit.generation
                result.total_costs += (
                    sum(result_unit.generation)
                    * result_unit.generator.marginal_costs
                )
            end
        end

        return result
    end
end

mutable struct ADMM
    iteration::Int
    gamma::Float64
    lambdas::Vector{Vector{Float64}}
    demand::Vector{Int}
    T::Vector{Int}
    generators::Vector{Generator}
    storages::Union{Vector{Storage}, Nothing}
    results::Vector{Result}
    converged::Bool

    function ADMM(gamma::Float64,
                  demand::Vector{Int},
                  generators::Vector{Generator},
                  storages::Union{Vector{Storage}, Nothing})
        admm = new()
        admm.iteration = 1
        admm.gamma = gamma
        admm.lambdas = [zeros(Float64, length(demand))]
        admm.demand = demand
        admm.T = collect(1:length(demand))
        admm.generators = generators
        admm.storages = storages
        admm.results = []
        admm.converged = false
        return admm
    end
end

function getSumOfIterationResults(iteration::Int)
    if length(admm.results) == 0
        empty_vector = zeros(Float64, length(demand))
        return empty_vector, empty_vector, empty_vector
    else
        return (
            admm.results[iteration].sum_G_t, 
            admm.results[iteration].sum_S_d_t,
            admm.results[iteration].sum_S_c_t
        )
    end
end

function getIterationResults(unit::Union{Generator, Storage}, iteration::Int)
    if length(admm.results) == 0
        empty_vector = zeros(Float64, length(demand))
        if typeof(unit) == Storage 
            return empty_vector, empty_vector
        else
            return empty_vector
        end
    else
        previous_results = admm.results[iteration]
        if typeof(unit) == Storage 
            return (
                previous_results.unit_to_result[unit].discharge,
                previous_results.unit_to_result[unit].charge
            )
        else
            return previous_results.unit_to_result[unit].generation
        end
    end
end
