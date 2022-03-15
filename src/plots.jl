export plotHistogram,
    plotHistogram_bus,
    plotHistogram_unc,
    plotHistogram_nodal,
    plotHistogram_branch,
    plotHistogram_6in9,
    plotHistogram_9in9

function plotHistogram(x::Vector; kwargs...)
    bins = haskey(kwargs, :bins) ? kwargs[:bins] : 70
    color = haskey(kwargs, :color) ? kwargs[:color] : "blue"
    alpha = haskey(kwargs, :alpha) ? kwargs[:alpha] : 0.5
    tight_layout() # spacing around labels
    xticks(rotation = 70) # rotate x axis labels
    plt[:hist](x, bins = bins, color = color, alpha = alpha)
    haskey(kwargs, :xlabel) ? xlabel(kwargs[:xlabel]) : nothing
    haskey(kwargs, :ylabel) ? ylabel(kwargs[:ylabel]) : nothing
    grid(true)
end

# Plots for generators or loads (2 buses)
function plotHistogram_bus(pg::Matrix, name::String, dir::String; kwargs...)
    haskey(kwargs, :fignum) ? figure(kwargs[:fignum]) : figure()
    plt[:subplot](221)
    plotHistogram(pg[1, :]; xlabel = name * "_1", ylabel = "ρ(" * name * "_1)", kwargs...)
    plt[:subplot](224)
    plotHistogram(pg[2, :]; xlabel = name * "_2", ylabel = "ρ(" * name * "_2)", kwargs...)
    println("Plotting: $dir/$name.pdf")
    savefig("$dir/$name.pdf")
    clf()
end

# Plot Histrogram for uncertainty bus parameters
function plotHistogram_unc(par1::Vector, par2::Vector, name::Vector{String}, dir::String; kwargs...)
    haskey(kwargs, :fignum) ? figure(kwargs[:fignum]) : figure()
    plt[:subplot](221)
    plotHistogram(par1; xlabel = name[1], ylabel = "ρ(" * name[1] * ")", kwargs...)
    plt[:subplot](224)
    plotHistogram(par2; xlabel = name[2], ylabel = "ρ(" * name[2] * ")", kwargs...)
    println("Plotting: $dir/$(name[1])_$(name[2]).pdf")
    savefig("$dir/$(name[1])_$(name[2]).pdf")
    clf()
end


# Plots for all nodes at once (4 buses)
function plotHistogram_nodal(v::Matrix, name::String, dir::String; kwargs...)
    haskey(kwargs, :fignum) ? figure(kwargs[:fignum]) : figure()
    plt[:subplot](221)
    plotHistogram(v[1, :]; xlabel = name * "_1", ylabel = "ρ(" * name * "_1)", kwargs...)
    plt[:subplot](222)
    plotHistogram(v[2, :]; xlabel = name * "_2", ylabel = "ρ(" * name * "_2)", kwargs...)
    plt[:subplot](224)
    plotHistogram(v[3, :]; xlabel = name * "_3", ylabel = "ρ(" * name * "_3)", kwargs...)
    plt[:subplot](223)
    plotHistogram(v[4, :]; xlabel = name * "_4", ylabel = "ρ(" * name * "_4)", kwargs...)
    println("Plotting: $dir/$name.pdf")
    savefig("$dir/$name.pdf")
    clf()
end

# Plots for all branches at once (5 branches)
function plotHistogram_branch(i::Matrix, name::String; kwargs...)
    haskey(kwargs, :fignum) ? figure(kwargs[:fignum]) : figure()
    plt[:subplot](332)
    plotHistogram(i[1, :]; xlabel = name * "_1", ylabel = "ρ(" * name * "_1)", kwargs...)
    plt[:subplot](335)
    plotHistogram(i[2, :]; xlabel = name * "_2", ylabel = "ρ(" * name * "_2)", kwargs...)
    plt[:subplot](334)
    plotHistogram(i[3, :]; xlabel = name * "_3", ylabel = "ρ(" * name * "_3)", kwargs...)
    plt[:subplot](336)
    plotHistogram(i[4, :]; xlabel = name * "_4", ylabel = "ρ(" * name * "_4)", kwargs...)
    plt[:subplot](338)
    plotHistogram(i[5, :]; xlabel = name * "_5", ylabel = "ρ(" * name * "_5)", kwargs...)
    println("Plotting: $dir/$name.pdf")
    savefig("plots/$name.pdf")
    clf()
end

# 6in 9x9 Histogram Plot
function plotHistogram_6in9(i::Matrix, name::String, dir::String; kwargs...)
    haskey(kwargs, :fignum) ? figure(kwargs[:fignum]) : figure()
    plt[:subplot](331)
    plotHistogram(i[1, :]; xlabel = name * "_1", ylabel = "ρ(" * name * "_1)", kwargs...)
    plt[:subplot](333)
    plotHistogram(i[2, :]; xlabel = name * "_2", ylabel = "ρ(" * name * "_2)", kwargs...)
    plt[:subplot](334)
    plotHistogram(i[3, :]; xlabel = name * "_3", ylabel = "ρ(" * name * "_3)", kwargs...)
    plt[:subplot](336)
    plotHistogram(i[4, :]; xlabel = name * "_4", ylabel = "ρ(" * name * "_4)", kwargs...)
    plt[:subplot](337)
    plotHistogram(i[5, :]; xlabel = name * "_5", ylabel = "ρ(" * name * "_5)", kwargs...)
    plt[:subplot](339)
    plotHistogram(i[6, :]; xlabel = name * "_6", ylabel = "ρ(" * name * "_6)", kwargs...)
    println("Plotting: $dir/$name.pdf")
    savefig("$dir/$name.pdf")
    clf()
end

# 9x9 Histogram Plot
function plotHistogram_9in9(m::Matrix, name::String, dir::String; kwargs...)
    haskey(kwargs, :fignum) ? figure(kwargs[:fignum]) : figure()
    for (i, row) in enumerate(eachrow(m))
        sub = 330 + i
        plt[:subplot](sub)
        plotHistogram(m[i, :]; xlabel = name * "_$i", ylabel = "ρ(" * name * "_$i)", kwargs...)
    end
    println("Plotting: $dir/$name.pdf")
    savefig("$dir/$name.pdf")
    clf()
end
