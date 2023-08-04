using DelimitedFiles

using Plots
using DataFrames
using LaTeXStrings
using KernelDensity
using Distributions
using Interpolations

import ColorSchemes

using Printf


"""
    get_tickslogscale(lims; skiplog=false)
Return a tuple (ticks, ticklabels) for the axis limit `lims`
where multiples of 10 are major ticks with label and minor ticks have no label
skiplog argument should be set to true if `lims` is already in log scale.
"""
function get_tickslogscale(lims::Tuple{T, T}; skiplog::Bool=false) where {T<:AbstractFloat}
    mags = if skiplog
        # if the limits are already in log scale
        floor.(lims)
    else
        floor.(log10.(lims))
    end
    rlims = if skiplog; 10 .^(lims) else lims end

    total_tickvalues = []
    total_ticknames = []

    rgs = range(mags..., step=1)
    for (i, m) in enumerate(rgs)
        if m >= 0
            tickvalues = range(Int(10^m), Int(10^(m+1)); step=Int(10^m))
            ticknames  = vcat([string(round(Int, 10^(m)))],
                              ["" for i in 2:9],
                              [string(round(Int, 10^(m+1)))])
        else
            tickvalues = range(10^m, 10^(m+1); step=10^m)
            ticknames  = vcat([string(10^(m))], ["" for i in 2:9], [string(10^(m+1))])
        end

        if i==1
            # lower bound
            indexlb = findlast(x->x<rlims[1], tickvalues)
            if isnothing(indexlb); indexlb=1 end
        else
            indexlb = 1
        end
        if i==length(rgs)
            # higher bound
            indexhb = findfirst(x->x>rlims[2], tickvalues)
            if isnothing(indexhb); indexhb=10 end
        else
            # do not take the last index if not the last magnitude
            indexhb = 9
        end

        total_tickvalues = vcat(total_tickvalues, tickvalues[indexlb:indexhb])
        total_ticknames = vcat(total_ticknames, ticknames[indexlb:indexhb])
    end

    return (total_tickvalues, total_ticknames)
end

"""
    fancylogscale!(p; forcex=false, forcey=false)
Transform the ticks to log scale for the axis with scale=:log10.
forcex and forcey can be set to true to force the transformation
if the variable is already expressed in log10 units.
"""
function fancylogscale!(p::Plots.Subplot; forcex::Bool=false, forcey::Bool=false)
    kwargs = Dict()
    for (ax, force, lims) in zip((:x, :y), (forcex, forcey), (xlims, ylims))
        axis = Symbol("$(ax)axis")
        ticks = Symbol("$(ax)ticks")

        if force || p.attr[axis][:scale] == :log10
            # Get limits of the plot and convert to Float
            ls = float.(lims(p))
            ts = if force
                (vals, labs) = get_tickslogscale(ls; skiplog=true)
                (log10.(vals), labs)
            else
                get_tickslogscale(ls)
            end
            kwargs[ticks] = ts
        end
    end

    if length(kwargs) > 0
        plot!(p; kwargs...)
    end
    p
end
fancylogscale!(p::Plots.Plot; kwargs...) = (fancylogscale!(p.subplots[1]; kwargs...); return p)
fancylogscale!(; kwargs...) = fancylogscale!(plot!(); kwargs...)

function mrbounddf(file ; header=false)
    dat = readdlm(file, header=header)
    df = DataFrame()
    df.r = Float64.(dat[:, 1])
    df.m = Float64.(dat[:, 2])

    return df
end

function lovembounddf(file ; header=false)
    dat = readdlm(file, header=header)
    df = DataFrame()
    df.m = Float64.(dat[:, 1])
    df.lambda = Float64.(dat[:, 2])

    return df
end

function getxy2(df)
    x = df.m
    y = df.lambda
    append!(x, df.m[1])
    append!(y, df.lambda[1])

    return (x, y)
end

function getxy(df)
    x = df.r
    y = df.m
    append!(x, df.r[1])
    append!(y, df.m[1])

    return (x, y)
end

@userplot MRExpPlot
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

@recipe function f(kc::MRExpPlot ; xlims = (7, 21), ylims = (0.5, 3.2), xticks = [8, 10, 12, 14, 16, 18, 20],
                   minorgrid = true, minorgridalpha = 0.15, minorgridlinewidth = 0.5, minorgridstyle = :dot, minorticks = 4,
                   xguide = L"$R$ $($km$)$", yguide = L"$M$ $($M$_{\odot})$")

    df1 = mrbounddf("GW170817_EOS_insensitive_90confidence_R1[km]_m1[Msun].dat")
    x1, y1 = getxy(df1)

    df2 = mrbounddf("GW170817_EOS_insensitive_90confidence_R2[km]_m2[Msun].dat")
    x2, y2 = getxy(df2)

    df3 = mrbounddf("GW170817_spectralEOS_wMaxm_90confidence_R1[km]_m1[Msun].dat")
    x3, y3 = getxy(df3)

    df4 = mrbounddf("GW170817_spectralEOS_wMaxm_90confidence_R2[km]_m2[Msun].dat")
    x4, y4 = getxy(df4)

    #FIXME: this one is a hack, I don't know how to get the 90% credibility bounds
    dat5, header5 = readdlm("J0030_2spot_RM.txt", ' ', '\n', header=true)

    df5 = DataFrame()
    df5.r = Float64.(dat5[1:3000, 1])
    df5.m = Float64.(dat5[1:3000, 2])

    x = vec(df5.r)
    y = vec(df5.m)

    k, kx, ky, ls = getkdecontour(x, y, ((-3.0, 3.0), (-3.0, 3.0)), 10)

    df6 = mrbounddf("J0030_Amstd_ST_PST_90_R[km]_m[Msun].dat")
    x6, y6 = getxy(df6)

    df7 = mrbounddf("J0740_IllinoisMaryland__NICER_XMM_relative__90_R[km]_m[Msun].dat")
    x7, y7 = getxy(df7)

    df8 = mrbounddf("J0740_Amstd_NICERxXMM___XPSI_STU_NSX_FIH__90_R[km]_m[Msun].dat")
    x8, y8 = getxy(df8)

    color1 = ColorSchemes.roma10.colors[1]
    color2 = ColorSchemes.bamako10.colors[7]
    color3 = ColorSchemes.tokyo10.colors[1]
    color4 = ColorSchemes.tokyo10.colors[8]

    color5 = ColorSchemes.acton10.colors[8]

    seriestype := :path
    fillalpha := 0.3

    legend --> false

    @series begin
        linecolor := color1
        fill := (0, color1)

        (x1, y1)
    end

    @series begin
        linecolor := color1
        fill := (0, color1)

        (x2, y2)
    end

    @series begin
        linecolor := color2
        fill := (0, color2)

        (x3, y3)
    end

    @series begin
        linecolor := color2
        fill := (0, color2)

        (x4, y4)
    end

    @series begin
        seriestype := :contour
        levels := [ls[1], last(ls)]
        fill := true

        color := [color3, color3]
        colorbar := false
        alpha := 0.3

        (collect(k.x), collect(k.y), k.density')
    end

    @series begin
        linecolor := color4
        fill := (0, color4)

        (x6, y6)
    end

    @series begin
        linecolor := color3
        fill := (0, color3)

        (x7, y7)
    end

    @series begin
        linecolor := color4
        fill := (0, color4)

        (x8, y8)
    end

    @series begin
        seriestype := :hline
        fill := 2.67
        linecolor := color5
        fillcolor := color5
        fillalpha := 0.3

        ([2.5])
    end

    @series begin
        seriestype := :hline
        linecolor := color5

        ([2.67])
    end

    @series begin
        seriestype := :hline
        fill := 3.1
        linecolor := color5
        fillcolor := color5
        fillalphaa := 0.3

        ([2.98])
    end

    @series begin
        seriestype := :hline
        linecolor := color5

        ([3.1])
    end
end

@userplot LoveMPlot
@recipe function f(kc::LoveMPlot ; xguide = L"$M$ $($M$_{\odot})$", yguide = L"$\Lambda$")

    df1 = lovembounddf("GW170817_EOS_insensitive_LC_90confidence_m1[Msun]_Love1.dat")
    x1, y1 = getxy2(df1)

    df2 = lovembounddf("GW170817_EOS_insensitive_LC_90confidence_m2[Msun]_Love2.dat")
    x2, y2 = getxy2(df2)

    df3 = lovembounddf("GW170817_spectralEOS_LC_wMaxm_90confidence_m1[Msun]_Love1.dat")
    x3, y3 = getxy2(df3)

    df4 = lovembounddf("GW170817_spectralEOS_LC_wMaxm_90confidence_m2[Msun]_Love2.dat")
    x4, y4 = getxy2(df4)

    color1 = ColorSchemes.roma10.colors[1]
    color2 = ColorSchemes.bamako10.colors[7]

    yscale := :log10

    xlim := (1.0, 3.5)
    ylim := (2, 2000)
    grid := false

    seriestype := :path
    fillalpha := 0.3

    legend --> false

    @series begin
        linecolor := color1
        fill := (0, color1)

        (x1, y1)
    end

    @series begin
        linecolor := color1
        fill := (0, color1)

        (x2, y2)
    end

    @series begin
        linecolor := color2
        fill := (0, color2)

        (x3, y3)
    end

    @series begin
        linecolor := color2
        fill := (0, color2)

        (x4, y4)
    end
end

function main()
    gr(size=(800, 600))

    #TODO: put annotations
    # pl = mrexpplot()
    # savefig(pl, "mrexp_paper.png")

    println(get_ticks([1, 2, 3, 4, 5, 6], 10))

    pl = lovemplot()
    fancylogscale!()
    savefig(pl, "lovemexp_paper.png")
end
