import Base: size, rand
using LinearAlgebra
using Combinatorics
using StatsBase


struct Graph{T}
    V::Set{T}
    E::Set{Set{T}}
    A::Matrix{Float64}
    D::Matrix{Float64}
    function Graph(V::Set{T}, E::Set{Set{T}}) where T <: Integer
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

num_vertices(G::Graph) = length(G.V)
num_edges(G::Graph) = length(G.E)
size(G::Graph) = (num_vertices(G), num_edges(G))

adjacency_matrix(G::Graph) = G.A
function adjacency_matrix(V::Set{T},E::Set{Set{T}}) where T
    n,m = length(V), length(E)
    A = zeros(Int,n,n)
    for edge in E
        i,j = edge
        A[i,j] = 1
        A[j,i] = 1
    end
    Symmetric(A)
end

degree_matrix(G::Graph) = G.D
function degree_matrix(V::Set{T},E::Set{Set{T}}) where T
    degrees = zeros(Int,length(V))
    for (i,j) in E
        degrees[i] += 1
        degrees[j] += 1
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
