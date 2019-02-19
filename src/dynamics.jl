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

function odefun_nonlinear(F::Formation)
    function dynamics(ẋ,x,p,t)
        p,n = size(F)
        s = NetworkState(x,p,n)
        ṡ = NetworkState(ẋ,p,n)
        nonlinear_formation_dynamics!(ṡ,s,F.G,F)
    end
end


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

function nonlinear_formation_dynamics!(ẋ::NetworkState, x::NetworkState, G::Graph{T}, F::Formation) where T
    p,n = size(x)
    Ẋ,X = ẋ.X, x.X
    Ẋ[:] .= 0
    for e in G.E
        i,j = edge(e)
        diff = norm(X[:,j] - X[:,i])
        Ẋ[:,i] += (diff^2 - F.d[i,j]^2)*(X[:,j] - X[:,i])
        if T <: AbstractUndirectedEdge
            Ẋ[:,j] += (diff^2 - F.d[i,j]^2 )*(X[:,i] - X[:,j])
        end
    end
end
