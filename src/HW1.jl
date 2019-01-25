using DifferentialEquations, Plots, Statistics, Random

include("graph.jl")

V = [1,2,3]
E = [[1,2],[2,3]]
G = Graph(V,E)

function linear_consensus(G::Graph,tf=10.,show_plot=false)
    L = graph_laplacian(G)

    function dynamics(ẋ,x,p,t)
        copyto!(ẋ, -L*x)
    end

    num_robots, m = size(G)
    x0 = float(rand(num_robots)*5)
    tspan = (0., tf)
    prob = ODEProblem(dynamics, x0, tspan)
    sol = solve(prob)

    if show_plot
        plot(sol,label=["x$i" for i = 1:num_robots],legend=:none,title="$num_robots robots, $m edges")
        x0_avg = mean(x0)
        plot!(collect(tspan),[x0_avg, x0_avg],color=:black,linestyle=:dash,label="average")
    end
end



G1 = random_graph(5,6)
is_connected(G1)

linear_consensus(G1,5,true)
G2 = random_graph(10,15)
is_connected(G2)
linear_consensus(G2,5,true)

G3 = random_graph(10,30)
is_connected(G3)
linear_consensus(G3,5,true)
