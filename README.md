# DecentralOPF.jl

This repository contains the source code for the implementation of a decentralized algorithm that solves an optimal power flow with multiple time periods and energy storage resources. The algorithm was an outcome of the master thesis _"Modeling Decentralized Electricity Markets - Solving Multi-Period Optimal Power Flow using Alternating Direction Method of Multipliers"_ written at the chair "Wirtschafts- und Infrastrukturpolitik" at TU Berlin and supervised by Dr. Richard Weinhold. 

## Abstract

Modern electricity systems have undergone an enormous transformation process in the last decades. One primary driver has been climate change and the involved actions to reduce carbon emissions. More and more renewable, non-dispatchable energy resources like photovoltaic and wind generators were integrated into almost every national electricity network, increasing the number of market participants and opposing challenges to the transmission grid operators to sustain a reliable electricity supply. On top of these changes, the availability of high-performance technology at very low costs enables new digital innovations to be on the forerun. One of those innovations is the trend toward decentralized systems. The most famous example is undoubtedly the cryptocurrency Bitcoin which provides an alternative to the centralized banking system and showcases a way to conduct transactions without an intermediary. This thesis investigates whether it is possible to decentralize an optimal power flow calculation that is a prevalent task of every transmission system operator. The optimal power flow considers multi-periods and the integration of energy storage resources. Based on the Alternating Direction Method of Multipliers and a review of current papers related to decentralized electricity markets, a decentralized algorithm is developed that solves an optimal power flow without a central entity knowing all sensitive information about the market participants. All computation is done by the market participants and is exchanged via an information network. The decentralized algorithm is applied to a three node case study system, and the obtained results are compared to a centralized optimal power flow. The comparison yields that the results are nearly identical except for minor differences in the per mille range. Some convergence problems were faced while implementing the mathematical formulations. They were removed by adapting the algorithm. Finally, a decentralized algorithm could be established and published as an open-source package to solve an optimal power flow with multi time periods and energy storage resources. The derivation and implementation of this algorithm are thoroughly documented in this thesis.

## Getting started

To run the decentralized model, you have to clone the repository:

```cmd
git clone https://github.com/rockstaedt/DecentralOPF.jl.git
```

Change into the directory and activate your Julia REPL.

```cmd
cd DecentralOPF.jl
julia
```

Go into the Julia package manager (`]`) and activate and initiate the environment to install all packages. 

```
(@v1.6) pkg> activate .
(DecentralOPF.jl) pkg> instantiate
```

Leave the package manager and run either model with the following commands:

```julia
julia> include("src/opf_central_reference.jl")
julia> include("src/opf_admm_decentral.jl")
```

Please be aware that you need a valid Gurobi license. An open-source solver like Clp was not tested yet. However, feel free
to replace the solver with a different one. The case study is defined in the file `src/cases/three_node.jl`.