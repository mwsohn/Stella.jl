"""
	anova(df::DataFrame, contvar::Symbol, groupvar::Symbol; pval=false)

Produces an one-way ANOVA table by default. `contvar' is a continous variable and `groupvar' is
the group variable. If `table` is set to false, you will have a DataFrame with ANOVA table values
reutrned. The p-value is in row 1, column 6 of the returned DataFrame.
"""
function anova(_df::DataFrame, dep::Symbol, cat::Symbol; table = true)

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
    pstr = pvalue < 0.0001 ? "< 0.0001" : @sprintf("%.5f", pvalue)

    outdf = DataFrame(
        Source=["Between", "Within", "Total"],
        DF=[dfbetween, dfwithin, dftotal],
        SS=[ssbetween, sswithin, sstotal],
        MS=[msbetween, mswithin, ""],
        F=[F, missing, missing],
        P=[pstr, missing, missing]
    )

	if pval == false
    	println("\nOne-Way Analysis of Variance: ", dep, " by ", cat)
    	pretty_table(outdf; formatters=(ft_nomissing, ft_printf("%.4f", [3, 4, 5])))
	else
		return pvalue
	end

end
