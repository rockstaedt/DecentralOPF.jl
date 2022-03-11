struct ResultStorage
    storage::Storage
    discharge::Vector{Float64}
    charge::Vector{Float64}
    level::Vector{Float64}
    penalty_term::PenaltyTerm
    U::Matrix{Float64}
    K::Matrix{Float64}
end

struct ResultGenerator
    generator::Generator
    generation::Vector{Float64}
    penalty_term::PenaltyTerm
    U::Matrix{Float64}
    K::Matrix{Float64}
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

mutable struct Result
    unit_to_result::Dict
    node_to_result::Dict{Node, ResultNode}
    generation::Vector{Float64}
    discharge::Vector{Float64}
    charge::Vector{Float64}
    penalty_term::PenaltyTerm
    avg_U::Matrix{Float64}
    avg_K::Matrix{Float64}
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
        result.avg_U = zeros(Float64, length(admm.L), length(admm.T))
        result.avg_K = zeros(Float64, length(admm.L), length(admm.T))

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

            result.avg_U += result_unit.U
            result.avg_K += result_unit.K

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
        result.avg_U = 1/(counter_subproblems) .* result.avg_U
        result.avg_K = 1/(counter_subproblems) .* result.avg_K

        result.line_utilization = admm.ptdf * result.injection

        return result
    end
end