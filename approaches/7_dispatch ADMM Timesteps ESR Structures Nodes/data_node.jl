using LinearAlgebra

mutable struct Node
    name::String
    demand::Vector{Int}
    slack::Bool
end

mutable struct Generator
    name::String
    marginal_costs::Int
    max_generation::Int
    plot_color::String
    node::Node
end

mutable struct Storage
    name::String
    marginal_costs::Int
    max_power::Int
    max_level::Int
    plot_color::String
    node::Node
end

mutable struct Line
    name::String
    from::Node
    to::Node
    max_capacity::Int
    susceptance::Int
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
    node_to_sum_penalty_t::Dict{Node, Vector{Float64}}
    total_costs::Float64
    sum_penalty_t::Vector{Float64}
    
    function Result(unit_to_result::Dict)
        result = new()

        result.unit_to_result = unit_to_result

        result.sum_G_t = zeros(Float64, length(admm.T))
        result.sum_S_d_t = zeros(Float64, length(admm.T))
        result.sum_S_c_t = zeros(Float64, length(admm.T))

        result.node_to_sum_penalty_t = Dict()
        for node in admm.nodes
            result.node_to_sum_penalty_t[node] = zeros(Float64, length(admm.T))
        end

        result.total_costs = 0
        result.sum_penalty_t = zeros(Float64, length(admm.T))

        for (unit, result_unit) in unit_to_result
            result.sum_penalty_t += result_unit.penalty_term
            result.node_to_sum_penalty_t[unit.node] += result_unit.penalty_term

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
    lambdas_n_t::Vector{Matrix{Float64}}
    T::Vector{Int}
    N::Vector{Int}
    nodes::Union{Vector{Node}, Nothing}
    generators::Vector{Generator}
    storages::Union{Vector{Storage}, Nothing}
    lines::Union{Vector{Line}, Nothing}
    results::Vector{Result}
    converged::Bool
    ptdf::Matrix{Float64}
    total_demand_t::Vector{Float64}
    node_to_demand::Dict{Node, Vector{Int}}
    node_to_id::Dict{Node, Int}
    node_to_units::Dict{Node, Vector{Union{Generator, Storage}}}

    function ADMM(gamma::Float64,
                  nodes::Union{Vector{Node}, Nothing},
                  generators::Vector{Generator},
                  storages::Union{Vector{Storage}, Nothing},
                  lines::Union{Vector{Line}, Nothing})
        admm = new()
        admm.iteration = 1
        admm.gamma = gamma
        admm.T = collect(1:length(nodes[1].demand))
        admm.N = collect(1:length(nodes))
        admm.lambdas_n_t = [zeros(Float64, length(admm.N), length(admm.T))]
        admm.nodes = nodes
        admm.generators = generators
        admm.storages = storages
        admm.lines = lines
        admm.results = []
        admm.converged = false
        admm.ptdf = calculate_ptdf(admm.nodes, admm.lines)
        admm.total_demand_t = zeros(length(admm.T))
        admm.node_to_demand = Dict()
        admm.node_to_id = Dict()
        admm.node_to_units = Dict()
        for (id, node) in enumerate(admm.nodes)
            admm.total_demand_t += node.demand
            admm.node_to_demand[node] = node.demand
            admm.node_to_id[node] = id
        end
        for unit in vcat(admm.generators, admm.storages)
            if haskey(admm.node_to_units, unit.node)
                push!(admm.node_to_units[unit.node], unit)
            else
                admm.node_to_units[unit.node] = [unit]
            end
        end
        return admm
    end
end

function getSumOfIterationResults(iteration::Int)
    if length(admm.results) == 0
        empty_vector = zeros(Float64, length(admm.T))
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
        empty_vector = zeros(Float64, length(admm.T))
        if typeof(unit) == Storage 
            return empty_vector, empty_vector
        else
            return empty_vector
        end
    else
        iteration_results = admm.results[iteration]
        if typeof(unit) == Storage 
            return (
                iteration_results.unit_to_result[unit].discharge,
                iteration_results.unit_to_result[unit].charge
            )
        else
            return iteration_results.unit_to_result[unit].generation
        end
    end
end

function calculate_ptdf(nodes::Vector{Node}, lines::Vector{Line})

    N = collect(1:length(nodes))
    SLACK = nothing

    incidence = zeros(Int, (length(lines), length(nodes)))
    susceptances = []

    for (l, line) in enumerate(lines)
        push!(susceptances, line.susceptance)
        for (n, node) in enumerate(nodes)
            if isnothing(SLACK)
                if node.slack
                    SLACK = n
                end
            end
            if line.from == node
                incidence[l, n] = 1
            elseif line.to == node
                incidence[l, n] = -1
            else
                incidence[l, n] = 0
            end
        end
    end
    
    B = Diagonal([susceptances...])
    # Line suceptance Matrix
    Bl = B*incidence
    # Nodes suceptance Matrix
    Bn = (incidence'*B)*incidence
    

    B_inv = zeros(length(N), length(N))

    B_inv[setdiff(N, SLACK), setdiff(N, SLACK)] = (
        inv(Bn[setdiff(N, SLACK), setdiff(N, SLACK)])
    )
    PTDF = Bl*B_inv
    return PTDF
end
