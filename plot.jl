using DelimitedFiles

using StatsPlots
# using Plots
gr(size=(800,600))
using DataFrames
using StatsBase
import KernelDensity
using Distributions
using Interpolations

function getkdecontour(x, y, clip, levels)
    m_x = median(x)
    m_y = median(y)

    dx_1 = m_x - quantile(x, 0.16)
    dx_h = quantile(x, 0.84) - m_x

    dy_1 = m_y - quantile(y, 0.16)
    dy_h = quantile(y, 0.84) - m_y

    xmin = m_x + clip[1][1] * dx_1
    xmax = m_x + clip[1][2] * dx_h

    ymin = m_y + clip[2][1] * dy_1
    ymax = m_y + clip[2][2] * dy_h

    k = KernelDensity.kde((x, y))
    kx = KernelDensity.kde(x)
    ky = KernelDensity.kde(y)

    ps = pdf.(Ref(k), x, y)

    ls = []
    for p in range(1.0 / levels, stop = 1 - 1.0 / levels, length = levels - 1)
        push!(ls, quantile(ps, p))
    end

    return (k, kx, ky, ls)
end

@userplot MRPlot

@recipe function f(kc::MRPlot; levels = 10, clip = ((-3.0, 3.0), (-3.0, 3.0)), color = ColorSchemes.roma10.colors[10],
                   xlimplot = (7, 15), ylimplot = (0.5, 3), xplotticks = [8, 10, 12, 14], yplotticks = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
                   minorticks = 4)
    x, y = kc.args
    x = vec(x)
    y = vec(y)

    k, kx, ky, ls = getkdecontour(x, y, clip, levels)

    legend --> false
    # layout := @layout [
    #     topdensity _
    #     contour{0.9w,0.9h} rightdensity
    # ]

    @series begin
        seriestype := :contour
        levels := [ls[1], last(ls)]
        fill := true
        # subplot := 2

        color := [color, color]
        colorbar := false
        alpha := 0.3

        xlims := xlimplot
        ylims := ylimplot
        xticks := xplotticks
        yticks := yplotticks

        minorgrid := true
        minorgridalpha := 0.15
        minorgridlinewidth := 0.5
        minorgridstyle := :dot
        minorticks := minorticks

        (collect(k.x), collect(k.y), k.density')
    end

    # ticks := nothing
    # xguide := ""
    # yguide := ""

    # @series begin
    #     seriestype := :density
    #     subplot := 1
    #     xlims := xlimplot
    #     ylims := (0, 1.1 * maximum(kx.density))
    #     color := color

    #     x
    # end

    # @series begin
    #     seriestype := :density
    #     subplot := 3
    #     orientation := :h
    #     xlims := (0, 1.1 * maximum(ky.density))
    #     ylims := ylimplot
    #     color := color

    #     y
    # end
end

function main()
    dat, header = readdlm("EoS-insensitive_posterior_samples.dat", ' ', '\n', header=true)

    df = DataFrame()
    df.m1 = Float64.(dat[:,1])
    df.m2 = Float64.(dat[:,2])
    df.lambda1 = Float64.(dat[:,3])
    df.lambda2 = Float64.(dat[:,4])
    df.r1 = Float64.(dat[:,5])
    df.r2 = Float64.(dat[:,6])

    color1 = ColorSchemes.roma10.colors[10]
    color2 = ColorSchemes.roma10.colors[1]

    pl = mrplot(df.r1, df.m1, color=color1)
    mrplot!(pl, df.r2, df.m2, color=color2)
    savefig(pl, "universalrel_mr_julia.png")
end
