import Base: size

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


function consensus_dynamics!(ẋ::NetworkState,x::NetworkState,G::Graph)
    p,n = size(x)
    X = x.X
    Ẋ = ẋ.X
    for i = 1:n
        Ni = neighbors(G,i)
        Ẋ[:,i] *= 0  # Make sure it is zero to start off with
        for j in neighbors(G,i)
            Ẋ[:,i] += X[:,j] - X[:,i]
        end
    end
end

function consensus_dynamics_laplacian!(ẋ::NetworkState,x::NetworkState,G::Graph)
    L = graph_laplacian(G)
    p,n = size(x)
    X = x.X
    Ẋ = ẋ.X
    for k = 1:p
        Ẋ[k,:] = -L*X[k,:]
    end
end


function simulate(N::Network, X0::NetworkState, dynamics::Function=odefun(N); tf::Real=10., show_plot=false)

    num_robots, m = size(N.G)
    tspan = (0., tf)
    prob = ODEProblem(dynamics, X0.x, tspan)
    sol = solve(prob)

    if show_plot
        plot(sol,label=["x$i" for i = 1:num_robots],legend=:none,title="$num_robots robots, $m edges")
        x0_avg = mean(x0)
        plot!(collect(tspan),[x0_avg, x0_avg],color=:black,linestyle=:dash,label="average")
    end
    return sol
end

get_trajectories(sol::DiffEqBase.AbstractODESolution,p,n) = [[[u[k,i] for u in [reshape(v,p,n) for v in sol.u]] for i = 1:n] for k = 1:p]

function plot(N::Network,sol; kwargs...)
    p,n = size(N)
    if p == 2
        X0 = reshape(sol.u[1],p,n)
        Xf = reshape(sol.u[end],p,n)
        X,Y = get_trajectories(sol,p,n)
        p = plot(X[1],Y[1],label="robot 1"; kwargs...)
        for i = 2:n
            plot!(X[i],Y[i],label="roboX0t $i")
        end
        scatter!(X0[1,:],X0[2,:],color=1:n,label="start",markershape=:square)
        scatter!(Xf[1,:],Xf[2,:],color=1:n,label="end",markershape=:circle)
    else
        p = plot(sol)
    end
    return p
end
