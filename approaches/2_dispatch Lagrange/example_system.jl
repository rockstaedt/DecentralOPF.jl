using DataFrames

include("utility.jl")

P = ["pv", "gas", "wind", "coal"]

plants = DataFrame(
    P = ["pv", "gas", "wind", "coal"],
    mc = [1, 50, 2, 30],
    gmax = [80, 150, 120, 300],
    node = ["N1", "N1", "N2", "N3"]
);

nodes = DataFrame(
    N = ["N1", "N2", "N3"],
    demand = [150, 30, 200],
    slack = [false, false, true]
);

lines = DataFrame(
    L = ["L1", "L2", "L3"],
    from = ["N2", "N3", "N2"],
    to = ["N1", "N1", "N3"],
    fmax = [40, 100, 100],
    susceptance = [1, 1, 2]

);

P = plants[:, :P]
N = nodes[:, :N]
L = lines[:, :L]

mc = Dict(k => v for (k,v) in eachrow(plants[:, [:P, :mc]]))
gmax = Dict(k => v for (k,v) in  eachrow(plants[:, [:P, :gmax]]))
map_np = Dict(k => v for (k,v) in eachrow(plants[:, [:P, :node]]))
map_pn = Dict(n => filter(row -> row.node == n, plants)[:, :P] for n in N)
demand =  Dict(k => v for (k,v) in eachrow(nodes[:, [:N, :demand]]))

ptdf = calculate_ptdf(nodes, lines)