using Plots, Statistics, Random, LinearAlgebra
include("graph.jl")
include("network.jl")
include("formation.jl")
include("dynamics.jl")
include("utils.jl")

# Simple 3-robot triangle
V = [1,2,3]
E = [[1,2],[2,3],[2,1],[2,3]]
G = Graph(V,E)
Xdes = Dict((1,2)=>[1.,1.],(2,3)=>[1.,-1.])
F = Formation(G,Xdes)
F.d
p,n = size(F)

X0 = [0. 0; 2 0.5; 0 2]'
s0 = NetworkState(X0)
sdot = NetworkState(zero(X0))
nonlinear_formation_dynamics!(sdot,s0,F.G,F)
dynfun = odefun_nonlinear(F)

sol = simulate(s0, dynfun, tf=10)
N = Network(G,dynfun,p,n)
sf = NetworkState(sol.u[end],N)

plot(N,sol)
plot!(sf, G)




# Create initials "BJ"
X_b = [0 0; 0 1; 0 2; 0 3; 0 4; 1 4; 2 4; 2.5 3; 2 2; 2.5 1; 2 0; 1 0; 1 2]
n_b = size(X_b,1)
X_j = [0 1; 0 0; 1 0; 2 0; 2 1; 2 2; 2 3; 2 4]
n_j = size(X_j,1)

# Add edges to formation
E = [[i,i+1] for i = 1:11]
push!(E,[12,1])  # Connect middle of "B"
push!(E,[3,13])  # Connect middle of "B"
push!(E,[13,9])  # Connect middle of "B"
E_j = [[i,i+1] for i = n_b.+(1:n_j-1)]
append!(E,E_j)
push!(E,[11,15])  # Connect letters at the bottom

# Create Graph
V = collect(1:n_b+n_j)   # Vertices
G = Graph(V,E)
p = 2                    # number of state dimensions
n,m = size(G)            # vertices, edges

# Create Xdes from positions and graph
bj = NetworkState([X_b; X_j .+ [3 0]]')
Xdes = get_edge_distances(bj,G)


# Create the "Formation" (location and graph)
F = Formation(G,Xdes)

# Set up problem and solve
s0 = NetworkState(rand(size(X)...))
dynfun = odefun_nonlinear(F)
sol = simulate(s0, dynfun, tf=100)
sf = NetworkState(sol.u[end],p,n)

# Plots
N = Network(G,dynfun,p,n)
plot(N,sol,legend=:none)
plot!(sf,G)

# Add edges
Random.seed!(3)
G = Graph(V,E)
add_random_edges!(G,40)
num_edges(G)
Xdes = get_edge_distances(bj,G)
F = Formation(G,Xdes)

s0 = NetworkState(rand(size(X)...))
dynfun = odefun_nonlinear(F)
sol = simulate(s0, dynfun, tf=20)
sf = NetworkState(sol.u[end],p,n)
maximum(calc_distance_error(sf,F))
plot(sf,aspect_ratio=:equal)
plot(sf,G)
