X = reshape(collect(1:12),3,4)
state = NetworkState(X)
state.x[1] = 10
@test state.X[1] == 10
state.X[1,2] = 12
@test state.x[4] == 12

state = NetworkState(2,4)
@test state.p == 2
@test state.n == 4
state.x[1] = 10
@test state.X[1] == 10


V = [1,2,3]
E = [[1,2],[2,3]]
G = Graph(V,E)

p = 1  # 2-D robots
n = num_vertices(G)
ẋ = NetworkState(p,n)
X0 = rand(1:10,p,n)
x = NetworkState(X0)
consensus_dynamics!(ẋ,x,G)
L = graph_laplacian(G)
@test -L*x.x == vec(ẋ.X )

p = 2  # 2-D robots
n = num_vertices(G)
ẋ = NetworkState(p,n)
X0 = rand(1:10,p,n)
x = NetworkState(X0)
N = Network(G,consensus_dynamics_laplacian!,p,n)

consensus_dynamics!(ẋ,x,G)
consensus_dynamics_laplacian!(ẋ2,x,G)
@test ẋ.X == ẋ2.X
