import Base: size

struct Formation
    G::Graph{T} where T
    Xdes::Dict{T,Vector{Float64}} where T
    Δ::Matrix{Float64}
    d::Matrix{Float64}
    function Formation(G::Graph{T2}, Xdes::Dict{T2,Vector{Float64}}) where T2 <: AbstractDirectedEdge
        if length(Xdes) != num_edges(G)
        error("Xdes must be the length of the edge set of G")
        end
        Δ,d = calcΔ(G,Xdes)
        new(G,Xdes,Δ,d)
    end
end

function Formation(G::Graph{T}, Xdes::Dict{NTuple{2,Int},Vector{Float64}}) where T <: AbstractUndirectedEdge
    E = G.E
    DE = Set{DirectedEdge}()
    Xdes2 = Dict{DirectedEdge,Vector{Float64}}()
    for e in E
        i,j = Tuple(edge(e))
        ij = (i,j)
        if (i,j) ∈ keys(Xdes)
            ij = (i,j)
            ji = (j,i)
        elseif (j,i) ∈ keys(Xdes)
            ij = (j,i)
            ji = (i,j)
        else
            error("Edge in graph doesn't have an associated vector")
        end

        dedge1 = DirectedEdge(ij)
        dedge2 = DirectedEdge(ji)
        push!(DE,dedge1)
        push!(DE,dedge2)
        Xdes2[dedge1] = Xdes[ij]
        Xdes2[dedge2] = -Xdes[ij]
    end
    G2 = Graph(G.V, DE)
    Formation(G2,Xdes2)
end

function size(F::Formation)
    e,xij = rand(Xdes)
    p = length(xij)
    n = num_vertices(F.G)
    return p,n
end

function calcΔ(G::Graph{T},Xdes) where T <: AbstractDirectedEdge
    n,m = size(G)
    e,xij = rand(Xdes)
    p = length(xij)
    Δ = zeros(p,n)
    d = zeros(n,n)
    for e in G.E
        i,j = edge(e)
        xij = Xdes[e]
        Δ[:,i] += xij
        d[i,j] = norm(Xdes[e])
    end
    return Δ,d
end

function calc_distance_error(x::NetworkState, F::Formation)
    diff = get_edge_distances(x, F.G)
    err = zeros(length(diff))
    for (k,e) in enumerate(F.G.E)
        i,j = edge(e)
        err[k] = abs(norm(diff[(i,j)]) - F.d[i,j])
    end
    return err
end

function get_edge_distances(x::NetworkState, G::Graph)
    p,n = size(x)
    X = x.X
    diff = zeros(num_edges(G),p)
    diff = Dict{NTuple{2,Int},Vector{Float64}}()
    for (k,e) in enumerate(G.E)
        i,j = edge(e)
        diff[(i,j)] = X[:,j] - X[:,i]
    end
    return diff
end
