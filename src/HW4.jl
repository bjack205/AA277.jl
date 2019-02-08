using Plots, Statistics, Random, LinearAlgebra
include("graph.jl")
include("formation.jl")
include("network.jl")
include("utils.jl")

# V = [1,2,3]
# E = [[1,2],[2,3],[2,1],[2,3]]
# A = [0 1 0; 1 0 1; 0 1 0]
# D = Diagonal([1,2,1])
# Edge.(E)
# G = Graph(V,E)
# DG = DGraph(V,E)
#
# size(G)
# size(DG)
# graph_laplacian(DG)
#
# V = Set(V)
# E = Set(DirectedEdge.(E))
# M = indicidence_matrix(V,E)
# D = degree_matrix(V,E,:in)
# A = adjacency_matrix(V,E)
# L = D-A
# Graph(V,E)
#
#
# e1 = WeightedDirectedEdge((1,2),2)
# e2 = WeightedDirectedEdge((2,3),3)
# e3 = WeightedDirectedEdge((2,1),4)
# e4 = WeightedDirectedEdge((2,1),5)
#
# V
# E = Set((e1,e2,e3,e4))
# DG = Graph(V,E)

V = [1,2,3]
E = [[1,2],[2,3],[2,1],[2,3]]
G = Graph(V,E)
Xdes = Dict((1,2)=>[1.,1.],(2,3)=>[1.,-1.])
F = Formation(G,Xdes)
p,n = size(F)

X0 = [0. 0; 0.5 1; 1 0]'
s0 = NetworkState(X0)
dynfun = odefun(F)
sol = simulate(s0, dynfun, tf=10)
N = Network(G,dynfun,p,n)
sf = NetworkState(sol.u[end],N)

plot(N,sol)
plot!(sf, G)
