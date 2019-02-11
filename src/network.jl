import Base: size
import Plots: plot, plot!
using DifferentialEquations

struct Network
    G::Graph    # Network graph
    f::Function  # Network Dynamics, f(ẋ::NetworkState, x::NetworkState, G::Graph)
    p::Int  # Number of states in each node
    n::Int  # Number of nodes
end
size(N::Network) = (N.p, N.n)

function odefun(N::Network)
    p,n = size(N)
    function dynamics(ẋ,x,p,t)
        s = NetworkState(x,N)
        ṡ = NetworkState(ẋ,N)
        N.f(ṡ,s,N.G)
    end
end

function proximity_consensus(p::Int,n::Int,R::Real,eps=1e-4)
    function dynamics(ẋ,x,params,t)
        s = NetworkState(x,p,n)
        ṡ = NetworkState(ẋ,p,n)
        V,E = proximity_graph(s,R,eps)
        L = Laplacian(V,E)
        consensus_dynamics_laplacian!(ṡ,s,L)
    end
end

"""
Generate ode function for formation control, where the graph is static and
    has the same edges as the target graph
"""
function odefun(F::Formation)
    function dynamics(ẋ,x,p,t)
        p,n = size(F)
        s = NetworkState(x,p,n)
        ṡ = NetworkState(ẋ,p,n)
        formation_dynamics!(ṡ,s,F.G,F)
    end
end



struct NetworkState{T}
    p::Int  # Number of states in each node
    n::Int  # Number of nodes
    x::AbstractVector{T}  # Vector of concatenated states (p*n,)
    X::AbstractMatrix{T}  # Matrix of states (p,n)
end

function NetworkState(X::AbstractMatrix{T}) where T
    p,n = size(X)
    x = vec(X)
    NetworkState{T}(p,n,x,X)
end

function NetworkState(p::Int,n::Int)
    X = zeros(p,n)
    NetworkState(X)
end

function NetworkState(x::Vector{T}, p::Int, n::Int) where T
    X = reshape(x,p,n)
    NetworkState(p,n,x,X)
end

NetworkState(x::Vector, N::Network) where {T} = NetworkState(x, N.p, N.n)

size(s::NetworkState) = (s.p,s.n)


consensus_dynamics_laplacian!(ẋ::NetworkState,x::NetworkState,G::Graph) = consensus_dynamics_laplacian!(ẋ,x,graph_laplacian(G))
function consensus_dynamics_laplacian!(ẋ::NetworkState,x::NetworkState,L::Laplacian)
    p,n = size(x)
    X = x.X
    Ẋ = ẋ.X
    for k = 1:p
        Ẋ[k,:] = -L*X[k,:]
    end
end

function formation_dynamics!(ẋ::NetworkState, x::NetworkState, G::Graph, F::Formation)
    consensus_dynamics_laplacian!(ẋ,x,G)
    p,n = size(x)
    Ẋ = ẋ.X
    for k = 1:p
        Ẋ[k,:] -= F.Δ[k,:]
    end
end



function simulate(N::Network, X0::NetworkState, dynamics::Function=odefun(N); tf::Real=10., show_plot=false)

    num_robots, m = size(N.G)
    sol = simulate(X0, dynamics, tf=tf)

    if show_plot
        plot(sol,label=["x$i" for i = 1:num_robots],legend=:none,title="$num_robots robots, $m edges")
        x0_avg = mean(x0)
        plot!(collect(tspan),[x0_avg, x0_avg],color=:black,linestyle=:dash,label="average")
    end
    return sol
end

function simulate(X0::NetworkState, dynamics::Function; tf::Real=10)
    tspan = (0., tf)
    prob = ODEProblem(dynamics, X0.x, tspan)
    solve(prob)
end

get_trajectories(sol::DiffEqBase.AbstractODESolution,p,n) = [[[u[k,i] for u in [reshape(v,p,n) for v in sol.u]] for i = 1:n] for k = 1:p]


function proximity_graph(s::NetworkState,R::Real,eps::Float64=1e-4)
    p,n = size(s)
    X = s.X
    V = Set(1:n)
    W = Set{WeightedEdge}()
    for (i,j) in combinations(1:n,2)
        dist = norm(X[:,i] - X[:,j])
        if dist <= R
            weight = 1/(dist+eps)
            push!(W,WeightedEdge((i,j),weight))
        end
    end
    V,W
end

function plot(N::Network,sol; kwargs...)
    p,n = size(N)
    if p == 2
        X0 = reshape(sol.u[1],p,n)
        Xf = reshape(sol.u[end],p,n)
        X,Y = get_trajectories(sol,p,n)
        p = plot(X[1],Y[1],label="robot 1"; kwargs...)
        for i = 2:n
        end
        scatter!(X0[1,:],X0[2,:],color=1:n,label="start",markershape=:square)
        scatter!(Xf[1,:],Xf[2,:],color=1:n,label="end",markershape=:circle)
    else
        p = plot(sol)
    end
    return p
end

function plot(state::NetworkState, G::Graph; kwargs...)
    p = plot()
    plot!(state, G; kwargs...)
    return p
end
function plot!(state::NetworkState, G::Graph; kwargs...)
    X = state.X
    for e in G.E
        (i,j) = edge(e)
        plot!(X[1,[i,j]],X[2,[i,j]],color=:black,label="")
    end
    scatter!(X[1,:],X[2,:],color=:blue,markersize=5,label=""; kwargs...)
end
