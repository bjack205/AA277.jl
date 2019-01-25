function plot_vertical_lines!(p::Plots.Plot,x::Vector; kwargs...)
	ylim = collect(ylims(p))
	plot_vertical_lines!(x, ylim; kwargs...)
end

function plot_vertical_lines!(x,ylim=[-100,100]; kwargs...)
    ys = [ylim for val in x]
    xs = [[val; val] for val in x]
    plot!(xs,ys, linestyle=:dash, color=:black, label=""; kwargs...)
end
