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

mutable struct PenaltyTerm
    energy_balance::Vector{Float64}
    flow1::Vector{Float64}
    flow2::Vector{Float64}
end

function sum_up(summand1::PenaltyTerm, summand2::PenaltyTerm)
    summand1.energy_balance += summand2.energy_balance
    summand1.flow1 += summand2.flow1
    summand1.flow2 += summand2.flow2
    return summand1
end

function get_empty_penalty_term()
    return PenaltyTerm(
        zeros(Float64, length(admm.T)),
        zeros(Float64, length(admm.T)),
        zeros(Float64, length(admm.T))
    )
end

struct ResultStorage
    storage::Storage
    discharge::Vector{Float64}
    charge::Vector{Float64}
    level::Vector{Float64}
    penalty_term::PenaltyTerm
    R_ref::Matrix{Float64}
    R_cref::Matrix{Float64}
end

struct ResultGenerator
    generator::Generator
    generation::Vector{Float64}
    penalty_term::PenaltyTerm
    R_ref::Matrix{Float64}
    R_cref::Matrix{Float64}
end

mutable struct ResultNode
    node::Node
    generation::Vector{Float64}
    discharge::Vector{Float64}
    charge::Vector{Float64}
    penalty_term::PenaltyTerm

    function ResultNode(node::Node)
        node_result = new()
        node_result.node = node
        node_result.generation = zeros(Float64, length(admm.T))
        node_result.discharge = zeros(Float64, length(admm.T))
        node_result.charge = zeros(Float64, length(admm.T))
        node_result.penalty_term = get_empty_penalty_term()
        return node_result
    end
end

function update(node_result::ResultNode,
                unit_result::Union{ResultStorage, ResultGenerator})
    if typeof(unit_result) == ResultStorage
        node_result.discharge += unit_result.discharge
        node_result.charge += unit_result.charge
    elseif typeof(unit_result) == ResultGenerator
        node_result.generation += unit_result.generation
    end
    node_result.penalty_term = sum_up(
        node_result.penalty_term,
        unit_result.penalty_term
    )
    return node_result
end

mutable struct Result
    unit_to_result::Dict
    node_to_result::Dict{Node, ResultNode}
    generation::Vector{Float64}
    discharge::Vector{Float64}
    charge::Vector{Float64}
    penalty_term::PenaltyTerm
    avg_R_ref::Matrix{Float64}
    avg_R_cref::Matrix{Float64}
    total_costs::Float64
    injection::Matrix{Float64}
    line_utilization::Matrix{Float64}
    
    function Result(unit_to_result::Dict)
        result = new()

        result.unit_to_result = unit_to_result

        result.generation = zeros(Float64, length(admm.T))
        result.discharge = zeros(Float64, length(admm.T))
        result.charge = zeros(Float64, length(admm.T))
        result.injection = zeros(Float64, length(admm.N), length(admm.T))
        result.avg_R_ref = zeros(Float64, length(admm.L), length(admm.T))
        result.avg_R_cref = zeros(Float64, length(admm.L), length(admm.T))

        result.node_to_result = Dict()
        for (node_id, node) in enumerate(admm.nodes)
            result.injection[node_id, :] -= node.demand
            result.node_to_result[node] = ResultNode(node)
        end

        result.total_costs = 0

        result.penalty_term = get_empty_penalty_term()

        for (unit, result_unit) in unit_to_result
            result.penalty_term = sum_up(
                result.penalty_term,
                result_unit.penalty_term
            )

            result.node_to_result[unit.node] = update(
                result.node_to_result[unit.node],
                result_unit
            )

            result.avg_R_ref += result_unit.R_ref
            result.avg_R_cref += result_unit.R_cref

            node_id = admm.node_to_id[unit.node]

            if typeof(unit) == Storage
                result.discharge += result_unit.discharge
                result.injection[node_id, :] += result_unit.discharge

                result.charge += result_unit.charge
                result.injection[node_id, :] -= result_unit.charge

                result.total_costs += result_unit.storage.marginal_costs * (
                    sum(result_unit.discharge + result_unit.charge)
                )
            elseif typeof(unit) == Generator
                result.generation += result_unit.generation
                result.injection[node_id, :] += result_unit.generation

                result.total_costs += (
                    sum(result_unit.generation)
                    * result_unit.generator.marginal_costs
                )
            end
        end

        # Calculate average of slack variables
        counter_subproblems = length(admm.generators) + length(admm.storages)
        result.avg_R_ref = 1/(counter_subproblems) .* result.avg_R_ref
        result.avg_R_cref = 1/(counter_subproblems) .* result.avg_R_cref

        result.line_utilization = admm.ptdf * result.injection

        return result
    end
end

mutable struct ADMM
    iteration::Int
    gamma::Float64
    lambdas::Vector{Vector{Float64}}
    mues::Vector{Matrix{Float64}}
    rhos::Vector{Matrix{Float64}}
    T::Vector{Int}
    N::Vector{Int}
    L::Vector{Int}
    nodes::Union{Vector{Node}, Nothing}
    generators::Vector{Generator}
    storages::Union{Vector{Storage}, Nothing}
    lines::Union{Vector{Line}, Nothing}
    results::Vector{Result}
    converged::Bool
    ptdf::Matrix{Float64}
    total_demand_t::Vector{Float64}
    node_id_to_demand::Dict{Int, Vector{Int}}
    node_to_id::Dict{Node, Int}
    node_to_units::Dict{Node, Vector{Union{Generator, Storage}}}
    L_max::Vector{Float64} 

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
        admm.L = collect(1:length(lines))
        admm.lambdas = [zeros(Float64, length(admm.T))]
        admm.mues = [zeros(Float64, length(admm.L), length(admm.T))]
        admm.rhos = [zeros(Float64, length(admm.L), length(admm.T))]
        admm.nodes = nodes
        admm.generators = generators
        admm.storages = storages
        admm.lines = lines
        admm.L_max = [line.max_capacity for line in admm.lines]
        admm.results = []
        admm.converged = false
        admm.ptdf = calculate_ptdf(admm.nodes, admm.lines)
        admm.total_demand_t = zeros(length(admm.T))
        admm.node_id_to_demand = Dict()
        admm.node_to_id = Dict()
        admm.node_to_units = Dict()
        if !isnothing(admm.nodes)
            for (id, node) in enumerate(admm.nodes)
                admm.total_demand_t += node.demand
                admm.node_id_to_demand[id] = node.demand
                admm.node_to_id[node] = id
            end
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

function get_results(iteration::Int)
    if length(admm.results) == 0
        empty_vector = zeros(Float64, length(admm.T))
        return empty_vector, empty_vector, empty_vector
    else
        return (
            admm.results[iteration].generation, 
            admm.results[iteration].discharge,
            admm.results[iteration].charge
        )
    end
end

function get_unit_results(unit::Union{Generator, Storage}, iteration::Int)
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

function get_node_results(iteration::Int, node::Node)
    if length(admm.results) == 0
        empty_vector = zeros(Float64, length(admm.T))
        return empty_vector, empty_vector, empty_vector
    else
        return (
            admm.results[iteration].node_to_result[node].generation, 
            admm.results[iteration].node_to_result[node].discharge,
            admm.results[iteration].node_to_result[node].charge
        )
    end
end

function get_slack_results(iteration::Int)
    if length(admm.results) == 0
        empty_matrix = zeros(Float64, length(admm.L), length(admm.T))
        return empty_matrix, empty_matrix
    else
        return (
            admm.results[iteration].avg_R_ref, 
            admm.results[iteration].avg_R_cref,
        )
    end
end

function get_node_id_to_result(iteration::Int)
    node_id_to_result = Dict()
    if length(admm.results) == 0
        for node in admm.nodes
            node_id_to_result[admm.node_to_id[node]] = ResultNode(node)
        end
    else
        result = admm.results[iteration]
        for (node, node_result) in result.node_to_result
            node_id_to_result[admm.node_to_id[node]] = node_result
        end
    end
    return node_id_to_result
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
