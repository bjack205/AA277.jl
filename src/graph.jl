import Base: size, rand, ==, convert, hash, length
import Plots.plot
using LinearAlgebra
using Combinatorics
using StatsBase
using Plots

AbstractUndirectedEdge  = AbstractSet{Int}
Edge = Set{Int}

struct WeightedEdge <: AbstractUndirectedEdge
    edge::Edge
    weight::Float64
    function WeightedEdge(e::Edge,w::Real)
        new(e::Edge,w)
    end
end
weight(e::WeightedEdge) = e.weight
weight(e::Edge) = 1
edge(e::WeightedEdge) = e.edge
edge(e::Edge) = e
WeightedEdge(e::NTuple{2,T},w::Real) where {T} = WeightedEdge(Edge(e),w)
(==)(e1::WeightedEdge,e2::WeightedEdge) = e1.edge == e2.edge
hash(e::WeightedEdge) = hash(e.edge)

abstract type AbstractDirectedEdge end
struct DirectedEdge <: AbstractDirectedEdge
    edge::NTuple{2,Int}
end
function DirectedEdge(e::AbstractVector{Int})
    if length(e) == 2
        DirectedEdge(Tuple(e))
    else
        error("Edge must have length of 2")
    end
end

struct WeightedDirectedEdge <: AbstractDirectedEdge
    edge::DirectedEdge
    weight::Float64
    function WeightedDirectedEdge(e::DirectedEdge, w::Real)
        new(e,w)
    end
end
WeightedDirectedEdge(e::NTuple{2,Int},w::Real) = WeightedDirectedEdge(DirectedEdge(e), w)
weight(e::WeightedDirectedEdge) = e.weight
weight(e::DirectedEdge) = 1
edge(e::WeightedDirectedEdge) = e.edge.edge
edge(e::DirectedEdge) = e.edge

AbstractEdge = Union{AbstractUndirectedEdge, AbstractDirectedEdge}
length(e::AbstractEdge) = length(edge(e))


struct Graph{T}
    V::Set{Int}
    E::Set{T}
    A::Matrix{Float64}
    D::Matrix{Float64}
    function Graph(V::Set{Int}, E::Set{T}) where T <: AbstractEdge
        for edge in E
            if length(edge) != 2
                delete!(E,edge)
            end
        end
        A = adjacency_matrix(V,E)
        D = degree_matrix(V,E)
        new{T}(V,E,A,D)
    end
end

function Graph(V::AbstractVector{T}, E::AbstractVector{TE}) where {T, TE <: AbstractVector}
    V = Set(V)
    E = Set(Edge.(E))
    Graph(V,E)
end

function DGraph(V::AbstractVector{T}, E::AbstractVector{TE}) where {T, TE <: AbstractVector}
    V = Set(V)
    E = Set(DirectedEdge.(E))
    Graph(V,E)
end



Laplacian{T} = Matrix{T}
# struct Laplacian{T<:Real}
#     A::Matrix{T}
# end

function Laplacian(V::Set{Int}, E::Set{T}) where {T<:AbstractEdge}
    Laplacian(degree_matrix(V,E) - adjacency_matrix(V,E))
end


num_vertices(G::Graph) = length(G.V)
num_vertices(L::Laplacian) = size(L,1)
num_edges(G::Graph) = length(G.E)
size(G::Graph) = (num_vertices(G), num_edges(G))

adjacency_matrix(G::Graph) = G.A
function adjacency_matrix(V::Set{Int},E::Set{T}) where {T<:AbstractEdge}
    n,m = length(V), length(E)
    A = zeros(n,n)
    for e in E
        i,j = edge(e)
        A[i,j] = weight(e)
        if T <: AbstractUndirectedEdge
            A[j,i] = weight(e)
        end
    end
    if T <: AbstractUndirectedEdge
        return Symmetric(A)
    else
        return A
    end
end

degree_matrix(G::Graph) = G.D
function degree_matrix(V::Set{Int},E::Set{T}) where {T<:AbstractEdge}
    degrees = zeros(length(V))
    for e in E
        i,j = edge(e)
        degrees[i] += weight(e)
        degrees[j] += weight(e)
    end
    Diagonal(degrees)
end
function degree_matrix(V::Set{Int},E::Set{T},direction::Symbol=:in) where {T<:AbstractDirectedEdge}
    degrees = zeros(length(V))
    for e in E
        i,j = edge(e)
        if direction == :in
            degrees[j] += weight(e)
        else
            degrees[i] += weight(e)
        end
    end
    Diagonal(degrees)
end

function indicidence_matrix(V::Set{Int}, E::Set{T}) where {T<:AbstractDirectedEdge}
    n,m = length(V), length(E)
    M = zeros(m,n)
    for (k,e) in enumerate(E)
        i,j = edge(e)
        w = weight(e)
        M[k,i] = w
        M[k,j] = -w
    end
    return M
end

graph_laplacian(G::Graph) = degree_matrix(G) - adjacency_matrix(G)
is_connected(G::Graph) = rank(graph_laplacian(G)) == num_vertices(G) - 1
random_graph(n::Int) = random_graph(n,rand(1:n^2))

function random_graph(n::Int,m::Int)
    V = 1:n
    E = collect(combinations(1:n,2))
    e_inds = sample(1:length(E),m,replace=false)
    E = E[e_inds]
    Graph(V,E)
end

function random_graph(vertices::AbstractVector=1:10,edges::AbstractVector=1:length(vertices)^2)
    n = rand(vertices)
    m = rand(edges)
    random_graph(n,m)
end

function neighbors(G::Graph, i::Int)
    A = adjacency_matrix(G)
    return findall(A[i,:].==1)
end



function Graph(L::Laplacian{T}) where {T}
    n = size(L,1)
    V = Set(1:n)
    W = Set{WeightedEdge}()
    for (i,j) in combinations(1:n,2)
        if L[i,j] < 0
            push!(W,WeightedEdge((i,j),-L[i,j]))
        end
    end
    Graph(V,W)
end
