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
