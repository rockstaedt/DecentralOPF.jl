using DataFrames

include("utility.jl")

P = ["pv", "gas", "wind", "coal"]

S = ["esr", "ps"]

plants = DataFrame(
    P = ["pv", "gas", "wind", "coal"],
    mc = [3, 50, 4, 30],
    gmax = [80, 150, 120, 300],
    node = ["N1", "N1", "N2", "N3"]
);

storages = DataFrame(
    P = ["esr", "ps"],
    mc = [1, 2],
    gmax = [10, 15],
    capa = [20, 30],
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
mc["esr"] = storages[1, :mc]
mc["ps"] = storages[2, :mc]
gmax = Dict(k => v for (k,v) in  eachrow(plants[:, [:P, :gmax]]))
gmax["esr"] = storages[1, :gmax]
gmax["ps"] = storages[2, :gmax]
capa = Dict("esr" => storages[1, :capa], "ps" => storages[2, :capa])
map_np = Dict(k => v for (k,v) in eachrow(plants[:, [:P, :node]]))
map_pn = Dict(n => filter(row -> row.node == n, plants)[:, :P] for n in N)
demand =  Dict(k => v for (k,v) in eachrow(nodes[:, [:N, :demand]]))

ptdf = calculate_ptdf(nodes, lines)