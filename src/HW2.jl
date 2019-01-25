using DifferentialEquations, Plots, Statistics, Random
import Plots: plot
include("graph.jl")
include("network.jl")
include("utils.jl")

Random.seed!(1)

function gen_switching_dynamics(Gs::Vector{Graph{TG}},T=1) where TG
    k = 1
    t1 = 0
    t0 = 0
    function switching_dynamics(ẋ,x,p,t)
        t == 0 ? t0 = 0 : nothing
        if t-t0 > T
            t0 = t
            k += 1
            k > length(Gs) ? k = 1 : nothing
        end
        s = NetworkState(x,N)
        ṡ = NetworkState(ẋ,N)
        N.f(ṡ,s,Gs[k])
    end
end

n = 5
X0 = float.(rand(1:n,p,n))
s0 = NetworkState(X0)
G1 = random_graph(n,6)
G2 = random_graph(n,4)
G3 = random_graph(n,5)
is_connected(G1)
is_connected(G2)
is_connected(G3)
Gs = [G1,G2,G3]
N = Network(G1,consensus_dynamics_laplacian!,p,n)

# Not switching
tf = 10
sol1 = simulate(N,s0,tf=tf)
plot(N,sol1,legend=:bottomleft,title="No Switching")

# Switching
T = 0.5
ode_switch = gen_switching_dynamics(Gs,T)
sol = simulate(N,s0,ode_switch,tf=tf)
plot(N,sol,legend=:bottomleft,title="Switching")

plt = plot(sol,title="Switching times")
plot_vertical_lines!(plt,collect(0:T:tf))
