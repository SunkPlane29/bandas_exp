using DelimitedFiles

using StatsPlots
# using Plots
gr(size=(800,600))
using DataFrames
using StatsBase
import KernelDensity
using Distributions
using Interpolations
using LaTeXStrings

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
                   xlimplot = (7, 21), ylimplot = (0.5, 3), xplotticks = [8, 10, 12, 14, 16, 18, 20],
                   yplotticks = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0], minorticks = 4)
    x, y = kc.args
    x = vec(x)
    y = vec(y)

    k, kx, ky, ls = getkdecontour(x, y, clip, levels)

    legend --> false

    @series begin
        seriestype := :contour
        levels := [ls[1], last(ls)]
        fill := true

        color := [color, color]
        colorbar := false
        alpha := 0.3

        xguide := L"$R$ $($km$)$"
        yguide := L"$M$ $($M$_{\odot})$"

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
end

import ColorSchemes

#TODO: export all this to hd5
function main()
    dat1, header1 = readdlm("EoS-insensitive_posterior_samples.dat", ' ', '\n', header=true)

    df1 = DataFrame()
    df1.m1 = Float64.(dat1[:,1])
    df1.m2 = Float64.(dat1[:,2])
    df1.lambda1 = Float64.(dat1[:,3])
    df1.lambda2 = Float64.(dat1[:,4])
    df1.r1 = Float64.(dat1[:,5])
    df1.r2 = Float64.(dat1[:,6])

    dat2, header2 = readdlm("Parametrized-EoS_maxmass_posterior_samples.dat", ' ', '\n', header=true)

    df2 = DataFrame()
    df2.m1 = Float64.(dat2[:, 1])
    df2.m2 = Float64.(dat2[:,2])
    df2.lambda1 = Float64.(dat2[:,3])
    df2.lambda2 = Float64.(dat2[:,4])
    df2.r1 = Float64.(dat2[:,5])
    df2.r2 = Float64.(dat2[:,6])

    # There is 40000 points in this posterior sample, it is expensive already with 3000 points
    dat3, header3 = readdlm("J0030_2spot_RM.txt", ' ', '\n', header=true)

    df3 = DataFrame()
    df3.r = Float64.(dat3[1:3000, 1])
    df3.m = Float64.(dat3[1:3000, 2])

    dat4 = readdlm("A_NICER_VIEW_OF_PSR_J0030p0451/ST_CST/run1/run1_nlive1000_eff0.3_noCONST_noMM_IS_tol-1.txt")

    df4 = DataFrame()
    df4.m = Float64.(dat4[1:10000, 4])
    df4.r = Float64.(dat4[1:10000, 5])

    dat5, header5 = readdlm("NICER+XMM-relative_J0740_RM.txt", ' ', '\n', header=true)

    df5 = DataFrame()
    df5.r = Float64.(dat5[1:3000, 1])
    df5.m = Float64.(dat5[1:3000, 2])

    color1 = ColorSchemes.roma10.colors[1]
    color2 = ColorSchemes.bamako10.colors[7]
    color3 = ColorSchemes.tokyo10.colors[1]
    color4 = ColorSchemes.tokyo10.colors[8]

    pl = mrplot(df1.r1, df1.m1, color=color1, levels=10)
    mrplot!(pl, df1.r2, df1.m2, color=color1, levels=10)
    mrplot!(pl, df2.r1, df2.m1, color=color2, levels=10)
    mrplot!(pl, df2.r2, df2.m2, color=color2, levels=10)
    mrplot!(pl, df3.r, df3.m, color=color3)
    mrplot!(pl, df4.r, df4.m, color=color4, levels=10)
    mrplot!(pl, df5.r, df5.m, color=color3)
    annotate!(pl, (9, 1.75, text("uni rel", 10, :center, :center, color1)))
    annotate!(pl, (10, 0.875, text("spec EoS", 10, :center, :center, color2)))
    annotate!(pl, (16.5, 1.875, text("IL/MD", 10, :center, :center, color3)))
    annotate!(pl, (14.5, 1.125, text("AM", 10, :center, :center, color4)))
    savefig(pl, "expmr_julia.png")
end
