using Test

V = [1,2,3]
E = [[1,2],[2,3]]
A = [0 1 0; 1 0 1; 0 1 0]
D = Diagonal([1,2,1])

G = Graph(V,E)
@test num_vertices(G) == 3
@test num_edges(G) == 2
@test adjacency_matrix(G) == A
@test degree_matrix(G) == D
@test graph_laplacian(G) == D-A

E2 = [[1,2],[2,3],[1,1],[2,1]]
G2 = Graph(V,E2)
@test num_edges(G2) == 2


e1 = WeightedEdge((2,1),4)
e2 = WeightedEdge((1,2),5)
e3 = WeightedEdge((1,3),3)
e1 == e2
W = Set((e1,e2,e3))
@test length(W) == 2

V = Set([1,2,3])
E = Set(Set.([[1,2],[2,3]]))
A = [0 1 0; 1 0 1; 0 1 0]
D = Diagonal([1,2,1])

adjacency_matrix(V,E)
@test_nowarn Laplacian(V,E)
@test_nowarn Laplacian(V,W)

degree_matrix(V,E)
degree_matrix(V,W)
adjacency_matrix(V,E)
adjacency_matrix(V,W)

for (i,j) in combinations(1:4,2)
    @show i,j
end

V = Set{WeightedEdge}()
push!(V,e1)
push!(V,e3)
