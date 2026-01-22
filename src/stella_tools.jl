"""
    smallest(da::AbstractArray, n::Int) or smallest(df::DataFrame, varname::Symbol, n::Int)

List the `n` smallest values in `da` in a descending order. The default `n` is 5.
"""
function smallest(da::AbstractArray; n::Int = 5)
	return sort(collect(skipmissing(da)))[1:n]
end
smallest(df::DataFrame,varname::Symbol; n = 5) = smallest(df[!,varname], n = n)

"""
    largest(da::AbstractArray, n::Int) or largest(df::DataFrame, varname::Symbol, n::Int)

List the `n` largest values in `da` in an ascending order. The default `n` is 5.
"""
function largest(da::AbstractArray; n::Int = 5)
	return sort(collect(skipmissing(da)))[end-n+1:end]
end
largest(df::DataFrame,varname::Symbol; n = 5) = largest(df[!,varname], n = n)

"""
    univariate(da::AbstractArray) or univariate(df::DataFrame,varname::Symbol)

Produce univariate statistics from the `da` and return a DataFrame with two columns,
:Statistic and :Value. Computed statistics include `N_Total` (number of rows
in the input DataArray), `N_Missing` (number of rows with missing values), `N_Used` (number of non-missing rows),
`Sum` (sum), `Mean` (mean), `SD` (standard deviation), `Var` (variance), `Min` (minimum),
`P25` (25th percentile), `Median` (median), `P75` (75th percentile), `Max` (maximum),
`Skewness` (skewness), and `Kurtosis` (kurtosis).
"""
function univariate(da::AbstractVector)

    da2 = dropmissing(da)

	return DataFrame(
		Statistic = [Symbol("N Total"),
			Symbol("N Missing"),
			Symbol("N Used"),
			:Sum,
			:Mean,
			:SD,
			:Var,
			:Min,
			:P25,
			:Median,
			:P75,
			:Max,
			:Skewness,
			:Kurtosis],
		Value = [size(da,1),
			size(da,1) - size(da2,1),
			size(da2,1),
			sum(da2),
			mean(da2),
			sqrt(var(da2)),
			var(da2),
			minimum(da2),
			quantile(da2,0.25),
			quantile(da2,0.5),
			quantile(da2,0.75),
			maximum(da2),
			skewness(da2),
			kurtosis(da2)]
	)
end
univariate(df::AbstractDataFrame,var::Symbol) = univariate(df[!,var])


"""
    univ(vec::AbstractVector)
    univ(df::AbstractDataFrame,varname::Union{Symbol,String})

Produce univariate statistics using `vec` or `df[!,varname]` and return
a DataFrame with univariate statistics and percentile values.
Computed statistics include `N Total` (number of rows in the input DataArray), 
`N Missing` (number of rows with missing values), `N Used` (number of non-missing rows),
`Sum` (sum), `Mean` (mean), `SD` (standard deviation), `Var` (variance), `Min` (minimum),
`P25` (25th percentile), `Median` (median), `P75` (75th percentile), `Max` (maximum),
`Skewness` (skewness), and `Kurtosis` (kurtosis).
"""
function univ(v::AbstractVector; table=true)

    v2 = collect(skipmissing(v)) # convert(Vector{Float64},v[.!ismissing.(v)])
    if length(v2) < 2
        error("Input vector has fewer than 2 non-missing values.")
    end

    S = smallest(v2)
    L = largest(v2)

    output = [
        "N"           size(v,1)                    "Minimum"               minimum(v2);
        "N Missing"   (size(v,1) - size(v2,1))     "5%"                    quantile(v2,0.05);
        "N Used"      size(v2,1)                   "10%"                   quantile(v2,0.1)
        "Sum"         sum(v2)                      "25%"                   quantile(v2, 0.25);
        "Mean"        mean(v2)                     "Median"                quantile(v2,0.5);
        "SD"          sqrt(var(v2))                "75%"                   quantile(v2, 0.75)
        "Var"         var(v2)                      "90%"                   quantile(v2, 0.9)
        "Skewness"    skewness(v2)                 "95%"                   quantile(v2,0.95);
        "Kurtosis"    kurtosis(v2)                 "Maximum"               maximum(v2);
        "Smallest"    ""                           "Largest"               "";
        1             S[1]                         1                       L[1];
        2             S[2]                         2                       L[2];
        3             S[3]                         3                       L[3];
        4             S[4]                         4                       L[4];
        5             S[5]                         5                       L[5];
    ]

    #@show output
    if table
        pretty_table(output, header = ["Statistic","","Percentile",""],hlines=[0,1,16],vlines=:none)
        return nothing
    end
    return output
end
function univ(df::AbstractDataFrame,s::Union{Symbol,String}; table=true)
    univ(df[!,s]; table=table)
end

function getdictval(dt::Dict,val)
    return haskey(dt,val) ? dt[val] : val
end

"""
    tabstat(df::DataFrame,varname::Symbol,groupvar::Symbol;s::Vector{Function} = [N,mean,sd,minimum,p25,median,p75,maximum ], skipmissing=false)

Produce a DataFrame that contains summary statistics for each `groupvar` subgroup
of the `varname` column in the `df`. The following are computed: `N` (total non-missing rows, default),
`mean` (mean), `sd` (standard deviation), `min` (minimum), `p25` (25th percentile),
`median` (median), `p75` (75th percentile), and `max` (maximum). If `table` is set to `false`,
the output will be returned as a DataFrame.
"""
function tabstat(indf::AbstractDataFrame, 
    var1::Union{String,Symbol}, 
    groupvar::Union{String,Symbol};
    s = Function[mean,sd,minimum,p25,median,p75,maximum], 
    skipmissing = true, 
    table = true)

    if length(s) == 0
        error("No statistic functions were specified.")
    end

    # strip prepended "Stella." from the function list
    namevec = [Symbol(replace(string(x),r"(Stella|Statistics)\." => "")) for x in s]

    # missing values in var1
    indf = indf[ismissing.(indf[!,var1]) .== false, :]

    # grouped df
    gdf = groupby(indf, groupvar, skipmissing = skipmissing)

    # groupvar and N
    outdf = combine(gdf, nrow => :n)

    # stats
    for j = 1:length(namevec) 
        outdf[!,namevec[j]] = [s[j](x[!,var1]) for x in gdf] 
    end

    # identify only rows with non-zero counts
    nz = outdf.n .!= 0

    if table
        pretty_table(outdf[nz, 2:end], 
        row_labels = outdf[nz,groupvar],
        row_label_column_title = label(indf,groupvar),
		show_subheader = false,
		vlines=[1])
        
        # anov = anova(indf, var1, groupvar)
        # println("One-way ANOVA: F(",anov.df[1],", ",anov.df[3],") = ",@sprintf("%.5f",anov.F[1]),", ", anov.pvalue[1] < 0.00001 ? "P < 0.00001" : @sprintf("P = %.5f",anov.pvalue[1]))
    else
        return outdf[nz,:]
    end
end

#----------------------------------------------------------------------------
# stats by subdataframe - begin
#----------------------------------------------------------------------------

"""
    p5(x::AbstractVector)
    p10(x::AbstractVector)
    p25(x::AbstractVector)
    p50(x::AbstractVector)
    p75(x::AbstractVector)
    p90(x::AbstractVector)
    p95(x::AbstractVector)

Produce percentile values (e.g., p5() produces 5% percentile, p10 10% percentile, etc)
from `x` AbstractVector. If `x` is a DataArray, it must be stripped of NA's before used in these functions.
Other percentile functions not shown here can easily be created by, for example,
`p15(x) = quantile(x,.15)`.
"""
p5(x) = quantile(x,.5)
p10(x) = quantile(x,.1)
p25(x) = quantile(x,.25)
p50(x) = quantile(x,.5)
p75(x) = quantile(x,.75)
p90(x) = quantile(x,.9)
p95(x) = quantile(x,.95)

"""
    cv(x::AbstractVector)

Is an alias of `variation()`.
"""
cv(x) = variation(x)

"""
    se(x::AbstractVector)

Is an alias of `sem(x)`.
"""
se(x) = sem(x)

"""
    sd(x::AbstractVector)

Is an alias of `std(x)`.
"""
sd(x) = std(x)

"""
    N(x::AbstractVector)

Is an alias of `size(x,1)`.
"""
N(x) = size(x,1)

"""
    eform(glmout::StatsModels.RegressionModel)

Produces exponentiated estimates, their standard errors, and 95% confidence intervals.
If you want to see odds ratios (ORs) or incidence rate ratios (IRRs)
after a logsitic regression or a Poisson regression, use this function to obtain
regression output instead of `coeftable`.

"""
function eform(glmout::StatsModels.TableRegressionModel)

    # family and link function
    cox = false
    if isa(glmout.model,GeneralizedLinearModel)
        distrib = glmout.model.rr.d
        linkfun = GLM.Link(glmout.model.rr)
    elseif isa(glmout.model, CoxModel)
	    cox = true
    else
        error("GLM or Cox model is required.")
    end

	coeftable2 = coeftable(glmout)

	# 95% CI
	if cox
        cv = quantile(Normal(),0.975)
        push!(coeftable2.cols, exp.(coeftable2.cols[1] .- cv .* coeftable2.cols[2]))
        push!(coeftable2.cols, exp.(coeftable2.cols[1] .+ cv .* coeftable2.cols[2])) 
    else
        coeftable2.cols[5] = exp.(coeftable2.cols[5])
        coeftable2.cols[6] = exp.(coeftable2.cols[6])
    end

	# estimates
	coeftable2.cols[1] = exp.(coeftable2.cols[1])

	# standard errors
	coeftable2.cols[2] = coeftable2.cols[1] .* coeftable2.cols[2]


	# rename column1 to OR
	if cox
	   coeftable2.colnms[1] = "HR"
       push!(coeftable2.colnms, "95% LB")
       push!(coeftable2.colnms, "95% UB")
	else
		coeftable2.colnms[1] = coeflab(distrib,linkfun)
	end

    # parse the row names and change variable names and values
    for i in 2:length(coeftable2.rownms)
        # parse row name and split into a tuple (varname, value)
        if occursin(": ",coeftable2.rownms[i])
            (varname,value) = split(coeftable2.rownms[i],": ")
        else
            (varname,value) = (coeftable2.rownms[i],"")
        end

        # If value is 1 and value label is Yes, it is a binary variable
        # do not print
    #     if value == 1 && ismatch(r"^ *yes *$"i,vlabel)
    #         coeftable2.rownms[i] = varlabel
    #     elseif value != ""
    #         coeftable2.rownms[i] = string(varlabel, ": ", vlabel)
    #     else
    #         coeftable2.rownms[i] = varlabel
    #     end
    end
    
	return coeftable2
end

function coeflab(d,l)
    if (isa(d,Bernoulli) || isa(d,Binomial)) && isa(l,LogitLink)
        return "OR"
    elseif (isa(d,Binomial) || isa(d,NegativeBinomial)) && isa(l,LogLink)
        return "RR"
    elseif isa(d,Poisson) && isa(l,LogLink)
        return "IRR"
    else
        return "exp(Est)"
    end
end

"""
    rcoeftable(m::StatisticalMModel, eform = false, robust = false; vce = :hc1, adjust = true)

Produces `coeftable` with robust standard errors and confidence intervals. Use `vce` option
to set the heteroskedasticity-consistent (HC) standard error functions computed using the
CovarianceMatrices.jl package. Choose the function using the following symbols:

    :hc0 for HC0()
    :hc1 for HC1()
    :hc2 for HC2()
    :hc3 for HC3()
    :hc4 for HC4()
    :hc4m for HC4m()
    :hc5 for HC5()

Adjust option is set to `true` by default to obtain confidence intervals that are consistent with
Stata or R values. If you want the original HC standard errors, set it to `false`.
"""
function rcoeftable(m::StatisticalModel; vce=:hc1, adjust=true)

    coeftab = coeftable(m)
    est = coeftab.cols[1]

    # column name
    coeftab.colnms[2] = "Robust SE"

    # cannot specify level in coeftable. It produces an error.
    level = 0.95

    # vce functions available
    vcedict = Dict(:hc0 => HC0(), :hc1 => HC1(), :hc2 => HC2(),
        :hc3 => HC3(), :hc4 => HC4(), :hc4m => HC4m(), :hc5 => HC5())

    # robust standard errors
    if typeof(m.model) <: LinearModel
        # use t-statistic
        se = stderror(HC1(), m)
        coeftab.cols[2] = se
        tval = est ./ se # t-statistic
        coeftab.cols[coeftab.teststatcol] = tval
        coeftab.cols[4] = 2 .* ccdf.(TDist(nobs(ols) - dof(ols) + 1), abs.(tval)) # pvalues
        cv = quantile.(TDist(nobs(ols) - dof(ols) - 1), [(1 - level) / 2, (1 + level) / 2])
        coeftab.cols[5] = est .+ tval .* cv[1] # lower
        coeftab.cols[6] = est .+ tval .* cv[2] # upper
        return coeftab
    end

    # GLM models
    # use Z-statistic
    se = stderror(vcedict[vce], m, prewhite=false)
    if adjust
        se = se .* (nobs(m) - 1) / nobs(m)
    end
    coeftab.cols[2] = se
    zval = est ./ se # z-statistic
    coeftab.cols[coeftab.teststatcol] = zval
    coeftab.cols[4] = 2 .* ccdf.(Normal(), abs.(zval)) # pvalues
    cv = quantile.(Normal(), [(1 - level) / 2, (1 + level) / 2])
    coeftab.cols[5] = est .+ se .* cv[1] # lower
    coeftab.cols[6] = est .+ se .* cv[2] # upper
    return coeftab
end

"""
    rconfint(m::StatisticalModel; vce = :hc1, adjust = true)

Produces confidence intervals based on the robust standard errors. Use `vce` option
to set the heteroskedasticity-consistent (HC) standard error functions computed using the
CovarianceMatrices.jl package. Choose the function using the following symbols:

    :hc0 for HC0()
    :hc1 for HC1()
    :hc2 for HC2()
    :hc3 for HC3()
    :hc4 for HC4()
    :hc4m for HC4m()
    :hc5 for HC5()

Adjust option is set to `true` by default to obtain confidence intervals that are consistent with
Stata or R values. If you want the original HC standard errors, set it to `false`.
"""
function rconfint(m::StatisticalModel; adjust=true, vce=:hc1)
    coeft = rcoeftable(m, adjust=adjust, vce=vce)
    return hcat(coeft.cols[5], coeft.cols[6])
end


#--------------------------------------------------------------------------
# pairwise correlations
#--------------------------------------------------------------------------
"""
    interleave(args...)

Produces a matrix whose rows were interleaved from all `args` matrices. All
matrices should be of the same dimensions. Only two-dimensional matrices are allowed.
"""
function interleave(args::AbstractMatrix...)
    if ndims(args[1]) != 2
        error("Only 2-dimensional arrays are allowed.")
    end
    len = size(args[1],1)
    nargs = length(args)
    dd = zeros(Float64,len*nargs,len)

    for i in 1:len
        for (j, mat) in enumerate(args)
            dd[ i*nargs - (nargs - j), :] = mat[i,:]
        end
    end
    dd
end



"""
    pwcorr(df::DataFrame, args::Symbol...; digits = 4, nobs = true, pvalue = true)
    pwcorr(df::DataFrame, args::Vector{Symbol}; digits = 4, nobs = true, pvalue = true)

Produces a pairwise (n x n) correlation matrix that contain correlation coefficients
between all `n` columns. It also displays the number of observations used for computing
each correlation coefficient and a p-value for testing H₀: r = 0. The number of digits
after the decimal point can be specified by `digits` option (default = 4). You can output
the number of observations used in the calculation, provide an option `nobs = true`. P-values
for each correlation coefficients can be output with `pvalue = true`. If you want
to have just the correlation coefficents, set `nobs = false` and `pvalue = false`.
"""
function pwcorr(indf::DataFrame, args::Vector{Symbol}; digits=4, nobs=false, pvalue=false)

    a = indf[!, args]
    colnames = names(a)
    cols = length(colnames)
    M = Matrix{Union{Missing,Any}}(missing, cols, cols)

    for i = 1:cols
        var1 = colnames[i]
        for j = 1:i
            if i == j
                M[i, i] = tuple(1.0, missing, sum(completecases(a[!, [colnames[j]]])))
                continue
            end
            var2 = colnames[j]
            b = dropmissing(a[:, [var1, var2]])
            if nrow(b) == 0
                error("No usable common observations in ", colnames[i], " and ", colnames[j])
            end
            ret = HypothesisTests.CorrelationTest(b[:, var1], b[:, var2])
            pval = 2 * min(ccdf(TDist(ret.n - 2), ret.t), ccdf(TDist(ret.n - 2), -ret.t))
            M[i, j] = tuple(ret.r, pval, ret.n)
        end
    end

    return PWCOR(M, colnames, digits, nobs, pvalue)
end
pwcorr(a::DataFrame, args::Symbol...; digits=4, nobs=false, pvalue=false) = pwcorr(a, [args...]; digits=digits, nobs=nobs, pvalue=pvalue)
struct PWCOR
    M::Matrix
    colnames::Vector
    digits::Int8
    nobs::Bool
    pvalue::Bool
end
function format_matrix(c)
    fmt = Printf.Format("%.$(c.digits)f")
    M = copy(c.M)
    for i = 1:size(M, 1)
        for j = 1:i
            ismissing(M[i, j]) && continue
            (r, p, n) = M[i, j]
            if c.pvalue && c.nobs
                M[i, j] = string(Printf.format(fmt, r), "\n", ismissing(p) ? "" : Printf.format(fmt, p), "\n", n, "\n")
            elseif c.pvalue
                M[i, j] = string(Printf.format(fmt, r), "\n", ismissing(p) ? "" : Printf.format(fmt, p), "\n")
            elseif c.nobs
                M[i, j] = string(Printf.format(fmt, r), "\n", n, "\n")
            else
                M[i, j] = Printf.format(fmt, r)
            end
        end
    end
    return M
end
function Base.show(io::IO, c::PWCOR)
    pretty_table(io, format_matrix(c),
        linebreaks=true,
        crop=:none,
        header=c.colnames,
        row_labels=c.colnames,
        formatters=(v, i, j) -> ismissing(v) ? "" : v,
        max_num_of_rows=10,
        hlines=[0, 1],
        vlines=[1])
end


#--------------------------------------------------------------------------
# ranksum(), signrank(), signtest()
#--------------------------------------------------------------------------
#
# struct ranksum_return
#     df::DataFrame
#     varname::Symbol
#     t::Float64
#     U::Float64
#     s²::Float64
#     var_t::Float64
#     z::Float64
#     pvalue::Float64
#     porder::Float64
# end
#
# function ranksum(df::DataFrame,varname::Symbol; group::Symbol = nothing)
#     if group == nothing
#         error("`group` variable is required")
#     end
#
#     df2 = df[completecases(df[[varname,group]]),[varname,group]]
#     df2[:__ranking] = tiedrank(df[varname])
#
#     df3 = by(df2,group) do subdf
#         DataFrame(
#             obs = size(subdf,1),
#             ranksum = sum(subdf[:__ranking]),
#             expected = size(subdf,1)*(size(df2,1)+1) / 2
#         )
#     end
#
#     # test statistic
#     t = df3[1,:ranksum]
#     n1 = df3[1,:obs]
#     U = t - n1*(n1+1) / 2
#     meanrank = mean(df2[:__ranking])
#     s² = var(df2[:__ranking])
#     vart = df3[1,:obs]*df3[2,:obs]*s² / size(df2,1)
#     z = (t - df3[1,:expected]) / sqrt(vart)
#     pval = ccdf()
#     porder = U / (df3[1,:obs]*df3[2,:obs])
#     pval = 2*Distributions.cdf(Distributions.Normal(), z)
#     return ranksum_return(df3,varname,t,U,s²,vart,z,pval,porder)
# end
# #
# # function print(r::ranksum_return)
# #     group = names(r.df)[1]
# #     print(r.df,"\n\n")
# #     @printf("         z = %.3f\n",r.z)
# #     @printf("Prob > |z| = %.3f\n",r.pvalue)
# #     @printf("P{%s(%s == %s)} >  P{%s(%s == %s)} = %.3f\n",
# #         string(r.varname),
# #         string(group),
# #         string(r.df[1,1]),
# #         string(r.varname),
# #         string(group),
# #         string(r.df[2,1]),
# #         r.porder
# #         )
# # end
#
# struct signrank_return
#     df::DataFrame
#     t::Float64
#     var_t::Float64
#     z::Float64
#     p::Float64
# end
#
# function signrank(df::DataFrame,var1::Symbol,var2::Symbol)
#
#     # complete cases only
#     df2 = df[completecases(df[!,[var1,var2]]),[var1,var2]]
#
#     # difference
#     d = df[!,var1] .- df[!,var2]
#
#     # sign
#     s = sign.(d)
#
#     # rank the absolute difference
#     r = tiedrank(abs(d)).*s
#
#     # test statistic
#     t = sum(r)
#
#     # dataframe output
#     n = size(df2,1)
#     n0 = sum(s .== 0)
#     et = (n*(n + 1)/2 - sum(s .== 0))/2
#     df = DataFrame(
#         sign = ["positive","negative","zero"],
#         obs = [
#             sum(s .== 1),
#             sum(s .== -1),
#             n0
#         ],
#         sum_ranks = [
#             sum((s .== 1) .* r),
#             abs(sum((s .== -1) .* r)),
#             n0
#         ],
#         expected = [et,et,n0]
#     )
#
#     varadjt = sum(r.^2) / 4
#     varunadjt = (n*(n + 1)*(2n + 1))/24
#     Δvarzeroadj = -1 * n0*(n0 + 1)*(2*n0 + 1) / 24
#     Δvartiesadj = varadjt - varunadjt - Δvarzeroadj
#
#     # z-value
#     z = (df[1,:sum_ranks] - et) / sqrt(varadjt)
#
#     # p-value
#     pval = 2*Distributions.cdf(Distributions.Normal(), z)
#
#     return signrank_return(df, et, varadjt, z, pval)
# end
# #
# # function print(sr::signrank_return)
# #     print(sr.df,"\n\n")
# #     @printf("Adj variance = %.3f\n",sr.var_t)
# #     @printf("           z = %.3f\n",sr.z)
# #     @printf("  Prob > |z| = %.3f\n",sr.p)
# # end
#
#
# struct signtest_return
#     df::DataFrame
#     p_left::Float64
#     p_right::Float64
#     p_both::Float64
# end
#
# function signtest(df::DataFrame,var1::Symbol,var2::Symbol)
#
#     # complete cases only
#     df2 = df[completecases(df[[var1,var2]]),[var1,var2]]
#
#     # difference
#     d = df[var1] .- df[var2]
#
#     # sign
#     s = sign.(d)
#
#     # test statistic
#     nplus = sum(s .> 0)
#     n0 = sum(s .== 0)
#     n_nonzero = size(df2,1) - n0
#     et = n_nonzero / 2
#
#     # output dataframe
#     df3 = DataFrame(
#         sign = ["positive","negative","zero"],
#         observed = [
#             sum(s .== 1),
#             sum(s .== -1),
#             n0
#         ],
#         expected = [et,et,n0]
#     )
#
#     # p1 - Ha: median of var1 - var2 > 0
#     p1 = Distributions.cdf(Distributions.Binomial(n_nonzero,.5),df3[2,:observed])
#
#     # p2 - Ha: median of var1 - var2 < 0
#     p2 = Distributions.cdf(Distributions.Binomial(n_nonzero,.5),df3[1,:observed])
#
#     # p2 - Ha: median of var1 - var2 != 0
#     p3 = min(1.0, 2.0*Distributions.cdf(Distributions.Binomial(n_nonzero,.5),df3[1,:observed]))
#
#     return signtest_return(df3, p1, p2, p3)
# end
# #
# # function print(sr::signtest_return)
# #     print(sr.df,"\n\n")
# #     @printf("Pr(#positive ≥ %d) = %.3f\n",sr.df[1,:observed],sr.p_left)
# #     @printf("Pr(#negative ≥ %d) = %.3f\n",sr.df[2,:observed],sr.p_right)
# #     @printf("Pr(#positive ≥ %d or #negative ≥ %d) = %.3f\n",sr.df[1,:observed],sr.df[2,:observed],sr.p_both)
# # end
# #

function dir(str::String)
  if occursin("*",str) || occursin("?",str)
    printdir(glob(str))
  elseif str == "." || str == ".."
    printdir(glob(string(str,"/*")))
  elseif str == ""
    printdir(readdir())
  else
    printdir(readdir(str))
  end
end
dir() = dir("")

function printdir(vstr::Vector{String})

  for i=1:length(vstr)
    vstr[i] = isdir(vstr[i]) ? string(vstr[i],"/") : vstr[i]
  end
  maxlen = maximum(length.(vstr))

  for i=1:length(vstr)
    println(rpad(vstr[i],maxlen),"  ", humanReadableByteCountBin(stat(vstr[i]).size))
  end
end

units   = ["K",  "M",  "G",  "T",  "P"]
bin_div = [0x1p10, 0x1p20, 0x1p30, 0x1p40, 0x1p40]

function humanReadableByteCountBin(bytes::Int64)
    b = bytes == typemin(Int64) ? typemax(Int64) : abs(bytes)
    if b < 1024
        return @sprintf("%d B", b)
    end
    sf = 40
    for i in 1:5
        if b < 0xfffcccccccccccc >> sf
            if i == 5
                b >>= 10
            end
            return @sprintf("%.1f %sB", b / bin_div[i], units[i])
        end
        sf -= 10
    end
    b = (b >> 20) / 0x1p40
    return @sprintf("%.1f EB", b)
end
function humansize(bytes::Int64)
    return humanReadableByteCountbin(bytes)
end

"""
	juliadate(::AbstractVector)
    juliadate(::Real)

converts SAS or stata date values to Julia Dates values.
"""
function juliadate(sasdt::AbstractVector)
    return Dates.Date.(Dates.UTD.(convert.(Int64,sasdt) .+ Dates.value(Date(1960,1,1))))
end
function juliadate(sasdt::Real)
    return juliadate([Int64(sasdt)])[1]
end

"""
	sasdate(::AbstractVector)
    sasdate(::Date)

Converts Julia Date values to SAS or stata values.
"""
function sasdate(juliadt::AbstractVector)
    return Dates.value.(juliadt) .- Dates.value(Date(1960, 1, 1))
end
function sasdate(juliadt::Date)
    return sasdate([juliadt])[1]
end


