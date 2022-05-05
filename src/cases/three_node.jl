node1 = Node("N1", [10, 250], false)
node2 = Node("N2", [50, 70], false)
node3 = Node("N3", [120, 200], true)

nodes = [node1, node2, node3]

line1 = Line("L1", node2, node1, 20, 1)
line2 = Line("L2", node3, node1, 45, 1)
line3 = Line("L3", node2, node3, 70, 2)

lines = [line1, line2, line3]

pv = Generator("pv", 3, 80, "yellow", node1)
wind = Generator("wind", 4, 120, "lightblue", node2)
coal = Generator("coal", 30, 300, "brown", node3)
gas = Generator("gas", 50, 120, "grey", node1)

generators = [pv, wind, coal, gas]

battery = Storage("battery", 1, 10, 20, "purple", node1)

storages = [battery]