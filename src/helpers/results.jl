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