using DifferentialEquations, Plots, Statistics, Random
import Plots: plot
include("graph.jl")
include("network.jl")
include("utils.jl")


n = 6
p = 2
X = reshape(rand(p*n)*2,p,n)
R = 1.2
state = NetworkState(X)
V,E = proximity_graph(state,R)
G = Graph(V,E)
is_connected(G)
plot(state,G)

# Static Network
N = Network(G,consensus_dynamics_laplacian!,p,n)
sol0 = simulate(N,state,tf=5)
plot(N,sol0,legend=:topleft,title="Static Graph")

# Changing Network
dynamics = proximity_consensus(p,n,R,0.1)
sol = simulate(state, dynamics, tf=4)
plot(N,sol,legend=:topleft,title="Proximity Graph")
sf = NetworkState(sol.u[end],p,n)
avg = mean(state.X,dims=2)
scatter!([avg[1]],[avg[2]],label="Initial Average",color=:black)


conn = zeros(length(sol.u))
for i = 1:length(sol.u)
    s = NetworkState(sol.u[i],p,n)
    V,S = proximity_graph(s,R)
    G = Graph(V,S)
    conn[i] = is_connected(G)
end
break_ind = findfirst(conn .== 0)
s = NetworkState(sol.u[break_ind],p,n)
avg2 = mean(s.X[:,[1,3,4,5,6]],dims=2)
scatter!([avg2[1]],[avg2[2]],label="Break Average",color=:red)

plot(sol0,legend=:none,title="Static Graph")
plot(sol,legend=:none,title="Proximity Graph")
