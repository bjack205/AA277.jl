import Base: size, rand, ==, convert, hash, length
import Plots.plot
using LinearAlgebra
using Combinatorics
using StatsBase
using Plots


AbstractEdge = AbstractSet{Int}
Edge = Set{Int}

struct WeightedEdge <: AbstractEdge
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
length(e::WeightedEdge) = length(e.edge)
WeightedEdge(e::NTuple{2,T},w::Real) where {T} = WeightedEdge(Edge(e),w)
(==)(e1::WeightedEdge,e2::WeightedEdge) = e1.edge == e2.edge
hash(e::WeightedEdge) = hash(e.edge)


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
    E = Set(Set.(E))
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
        A[j,i] = weight(e)
    end
    Symmetric(A)
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
