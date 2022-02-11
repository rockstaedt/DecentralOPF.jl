function print_results(print_generators::Bool,
    print_storages::Bool,
    print_penalty::Bool,
    iteration::Int)
unit_to_result = admm.results[iteration].unit_to_result

for (unit, result) in unit_to_result
if typeof(unit) == Generator && print_generators
print("$(result.generator.name): ")
println(result.generation)
elseif typeof(unit) == Storage && print_storages
println("$(result.storage.name): ")
println("\tdischarge: $(result.discharge)")
println("\tcharge: $(result.charge)")
println("\tlevel: $(result.level)")
end
end
# TODO: Split between nodes and lines
if print_penalty
println("Penalty Sum per Node:")
node_to_result = admm.results[iteration].node_to_result
for (node, node_result) in node_to_result
println("\t$(node.name): $(node_result.penalty_term)")
end
println("Penalty Sum Total: $(admm.results[iteration].penalty_term)")
end
end

function print_duals(iteration::Int)
println(
"\nlambdas: $(admm.lambdas[iteration])"
* "\nmues: $(admm.mues[iteration])"
* "\nrhos: $(admm.rhos[iteration])\n"
)
end