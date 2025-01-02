struct AOV
    title::Vector{String}
    ss::Vector{Float64}
    df::Vector{Int64}
    ms::Vector{Float64}
    F::Float64
    pvalue::Float64
end

"""
	anova(::DataFrame, contvar::Symbol, groupvar::Symbol)
    anova(::StatsModels.TableRegressionModel)

Produces an one-way ANOVA table. `contvar' is a continous variable and `groupvar' is
the group variable. 
"""
function anova(_df::DataFrame, dep::Symbol, cat::Symbol)

    # establish data
    df = _df[completecases(_df[:, [dep, cat]]), [dep, cat]]

    gdf = groupby(df, cat)
    groups = []
    for subdf in gdf
        push!(groups, subdf[!, dep])
    end

    # grand mean
    μ = mean(vcat(groups...))

    # levels and number of groups
    lev = sort(collect(values(gdf.keymap)))
    k = length(lev)

    # group means and between group sum of squares
    groupmean = mean.(groups)
    groupvar = var.(groups)
    n = length.(groups)
    ssbetween = sum(n .* (groupmean .- μ) .^ 2)
    sswithin = sum((n .- 1) .* groupvar)

    sstotal = ssbetween + sswithin
    dfbetween = k - 1
    dfwithin = sum(n) - k
    dftotal = sum(n) - 1
    mswithin = sswithin / dfwithin
    msbetween = ssbetween / dfbetween
    F = msbetween / mswithin

    pvalue = Distributions.ccdf(Distributions.FDist(dfbetween, dfwithin), F)
    pstr = pvalue < 0.0001 ? "< 0.0001" : @sprintf("%.4f", pvalue)

    return AOV(
        ["Between", "Within", "Total"],
        [ssbetween, sswithin, sstotal],
        [dfbetween, dfwithin, dftotal],
        [msbetween, mswithin, sstotal / dftotal],
        F,
        pvalue
    )

end
function anova(glmmodel)
    tss = nulldeviance(glmmodel)
    rss = deviance(glmmodel)
    mss = tss - rss
    tdf = nobs(glmmodel) - 1
    mdf = GLM.dof(glmmodel) - 1
    rdf = tdf - mdf
    F = (mss / mdf) / (rss / rdf)
    return AOV(
        ["Model", "Residual", "Total"],
        [mss, rss, tss],
        [mdf, rdf, tdf],
        [(tss - rss) / mdf, rss / rdf, tss / tdf],
        F,
        Distributions.ccdf(Distributions.FDist(mdf, rdf), F)
    )
end
import Base.show
function Base.show(io::IO, a::AOV)

    pstr = a.pvalue < 0.0001 ? "< 0.0001" : @sprintf("%.4f", a.pvalue)
    println("\nAnalysis of Variance\n")
    pretty_table(DataFrame(
            Source=a.title, SS=a.ss, DF=a.df, MS=a.ms,
            F=[a.F, missing, missing],
            P=[pstr, missing, missing]
        );
        formatters=(ft_nomissing, ft_printf("%.3f", [2, 4, 5])),
        hlines=[1, 3],
        vlines=[1],
        show_subheader=false
    )

end
