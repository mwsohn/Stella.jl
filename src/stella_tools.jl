dropmissing(x) = collect(skipmissing(x))

"""
    smallest(da::AbstractArray, n::Int) or smallest(df::DataFrame, varname::Symbol, n::Int)

List the `n` smallest values in `da` in a descending order. The default `n` is 5.
"""
function smallest(da::AbstractArray; n::Int = 5)
	return sort(collect(skipmissing(da)))[1:n]
end
smallest(df::DataFrame,varname::Symbol; n = 5) = smallest(df[varname], n = n)

"""
    largest(da::AbstractArray, n::Int) or largest(df::DataFrame, varname::Symbol, n::Int)

List the `n` largest values in `da` in an ascending order. The default `n` is 5.
"""
function largest(da::AbstractArray; n::Int = 5)
	return sort(collect(skipmissing(da)))[end-n+1:end]
end
largest(df::DataFrame,varname::Symbol; n = 5) = largest(df[varname], n = n)

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
		Statistic = [:N_Total,
			:N_Missing,
			:N_Used,
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
univariate(df::DataFrame,var::Symbol) = univariate(df[var])

function strval(val::AbstractFloat)
  return @sprintf("%.2f",val)
end
function strval(val::AbstractFloat, decimals::Int)
    if decimals == 1
        return @sprintf("%.1f",val)
    elseif decimals == 2
        return @sprintf("%.2f",val)
    elseif decimals == 3
        return @sprintf("%.3f",val)
    elseif decimals == 4
        return @sprintf("%.4f",val)
    elseif decimals == 5
        return @sprintf("%.5f",val)
    elseif decimals == 6
        return @sprintf("%.6f",val)
    else
        return @sprintf("%f",val)
    end
end
function strval(val::Integer)
  return @sprintf("%.0d",val)
end

function getdictval(dt::Dict,val)
    return haskey(dt,val) ? dt[val] : val
end

"""
    tabstat(df::DataFrame,varname::Symbol,groupvar::Symbol;s::Vector{Function} = [N,mean,sd,minimum,p25,median,p75,maximum ],wt::Union{Nothing,Symbol}=nothing)

Produce a DataFrame that contains summary statistics for each `groupvar` subgroup
of the `varname` column in the `df`. The following are computed: `n` (total non-missing rows),
`mean` (mean), `sd` (standard deviation), `min` (minimum), `p25` (25th percentile),
`median` (median), `p75` (75th percentile), and `max` (maximum).
"""
function tabstat(indf::DataFrame, var1::Symbol, groupvar::Symbol; s::Vector{Function} = [N,mean,sd,minimum,p25,median,p75,maximum ],wt::Union{Nothing,Symbol}=nothing)

    if length(s) == 0
        error("No statistic functions were specified.")
    end

    # prepend Stella to the functions
    namevec = [Symbol(replace(string(x),"Stella.","")) for x in s]

    # number of levels in the groupvar
    gvnum = length(DataFrames.levels(indf[groupvar]))

    outdf = DataFrame()
    outdf[groupvar] = Vector{Union{Missing,eltype(groupvar)}}(gvnum)
    for j = 1:length(namevec)

        # for N, create a vector of integers
        if namevec[j] == :N
            outdf[namevec[j]] = Vector{Union{Missing,Int}}(gvnum)
        # for all other stats, create float64 vectors
        else
            outdf[namevec[j]] = Vector{Union{Missing,Float64}}(gvnum)
        end
    end

    for subdf in groupby(indf, groupvar)
        gidx = findfirst(lev,subdf[1,groupvar])
        da = collect(skipmissing(subdf[var1]))
        if length(da) == 0
            continue
        end
        outdf[gidx,groupvar] = subdf[1,groupvar]
        for j = 1:length(namevec)
            outdf[gidx,namevec[j]] = s[j](da)
        end
    end

    sort!(outdf,groupvar)

    return outdf
end

#----------------------------------------------------------------------------
# stats by subdataframe - begin
#----------------------------------------------------------------------------

"""
    pickone(df::DataFrame,groupvars::Array{Symbol,1})

Create a DataArray that identifies one record in a `groupvar` in the `df` DataFrame.
This function creates a variable similar to the Stata command
`egen byte pickone = tag(groupvar)`. `groupvars` can be an array of Symbols or a single
Symbol.
"""
function pickone(df::DataFrame,groupvars::Array{Symbol,1})
    df[:___obs___] = collect(1:size(df,1))
    done = zeros(Int8,size(df,1))
    for subdf in groupby(df, groupvars)
        done[subdf[1,:___obs___]]=1
    end
    delete!(df,:___obs___)
    return convert(Vector{Union{Missing,Int8}},done)
end
pickone(df::DataFrame,groupvar::Symbol) = pickone(df, [groupvar])

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
# substat - pass functions to produce the stats by subframe
# functions: mean, var, std, sum, count, kurtosis, skewness, median, percentiles, maximum, minimum
#            cov, se


"""
    substat(df::DataFrame, varname::Symbol, groupvars::Array{Symbol,1}, func::Function)

Produces a DataArray of the same length as the `df` that contains summary statistics of
the `varname` column computed by `func` for each subgroup stratified by `groupvar`.
Multiple columns can be used as `groupvar`'s if provided in an array of Symbols.
Any function that can take one numeric vector as input can be used. Percentile functions
such as those that take two arguments (e.g., `quantile(x,.25)` for a 25th percentile of x)
should first be converted to a function that takes one argument (e.g, `p25(x) = quantile(x,.25)`).
The following percentile functions are available: `p5`, `p10`, `p25`, 'p50', `p75`, and `p90`.
This function emulates the Stata's `egen functions` such as `egen testmean = mean(test), by(groupvar)`.

```
julia> df[:testmean] = substat(df,:test, :class, mean)
```
"""
function substat(df::DataFrame, varname::Symbol, groupvars::Vector{Symbol}, func::Function)
    if (eltype(df[varname]) <: Number) == false
        error("Only numeric variables are allowed.")
    end

    df2 = by(df,groupvars) do subdf
        da = Vector(collect(skipmissing(subdf[varname])))
        if length(da) == 0
            DataFrame(x1 = NA)
        else
            DataFrame(x1 = func(da))
        end
    end

    df4 = DataFrame()
    df4[groupvars] = df[groupvars]
    df4[:___obs___] = collect(1:size(df4,1))

    df5 = join(df4, df2, on = groupvars, kind = :left)
    sort!(df5,cols=[:___obs___])
    return df5[:x1]
end
substat(df::DataFrame, varname::Symbol, groupvar::Symbol, func::Function) = substat(df,varname,[groupvar],func)

#
# #----------------------------------------------------------------------------
# # stats by subdataframe - end
# #----------------------------------------------------------------------------
# """
#     recode(da::DataArray,dict::Dict; restna=false)
#     recode(df::DataFrame,v::Symbol,dict::Dict; restna=false)
#
# Recode values in `da` or `df[v]` to values specified in `dict` dictionary.
# An option `restna = true` will convert all values in the original data array that are
# not a key in `dict` to NAs. Values in `dict` must be the same type as the original
# value or integer.
#
# ## Example
#
# ```
# julia> df = DataFrame(race = ["White","White","Black","Other","Hispanic"], sex = ["M","F","M","M","F"])
# 5×2 DataFrames.DataFrame
# │ Row │ race       │ sex │
# ├─────┼────────────┼─────┤
# │ 1   │ "White"    │ "M" │
# │ 2   │ "White"    │ "F" │
# │ 3   │ "Black"    │ "M" │
# │ 4   │ "Other"    │ "M" │
# │ 5   │ "Hispanic" │ "F" │
#
# julia> df[:race2] = recode(df,:race,Dict("White" => 1,"Black" => 2, "Hispanic" => 3, "Other" => 4))
# 5-element DataArrays.DataArray{Int8,1}:
#  1
#  1
#  2
#  4
#  3
#
#  julia> df
# 5×3 DataFrames.DataFrame
# │ Row │ race       │ sex │ race2 │
# ├─────┼────────────┼─────┼───────┤
# │ 1   │ "White"    │ "M" │ 1     │
# │ 2   │ "White"    │ "F" │ 1     │
# │ 3   │ "Black"    │ "M" │ 2     │
# │ 4   │ "Other"    │ "M" │ 4     │
# │ 5   │ "Hispanic" │ "F" │ 3     │
#
# ```
#
# """
# function recode(da::AbstractArray, coding::Dict; restna = false)
#     val = values(coding)
#
#     # if the da is not integer type
#     # check to see if all values in the coding dictionary are integers or NAs
#     # if so, construct a return data array whose elements are integers
#     # otherwise, keep the original data type
#     if !(eltype(da) <: Integer) && sum([typeof(v) <: Integer || isna(v) for v in val]) == length(val)
#         ra = DataArray(Int64,length(da))
#     else
#         ra = DataArray(eltype(da),length(da))
#     end
#
#     for i in 1:length(da)
#         if isna(da[i])
#             continue
#         end
#         if restna
#             ra[i] = haskey(coding,da[i]) ? coding[da[i]] : NA
#         else
#             ra[i] = haskey(coding,da[i]) ? coding[da[i]] : da[i]
#         end
#     end
#     if eltype(ra) <: Integer
#         return compress(ra)
#     end
#     return ra
# end
# recode(df::DataFrame,varname::Symbol,coding::Dict; restna = false) = recode(df[varname],coding,restna=restna)

#----------------------------------------------------------------------------
# eform
#----------------------------------------------------------------------------
function eform(coeftbl::StatsBase.CoefTable)
	coeftable2 = coeftbl

	# estimates
	coeftable2.cols[1] = exp.(coeftable2.cols[1])

	# standard errors
	coeftable2.cols[2] = coeftable2.cols[1] .* coeftable2.cols[2]

	# rename column1 to OR
	coeftable2.colnms[1] = "OR"

	return coeftable2
end

function eform(coeftbl::StatsBase.CoefTable, labels::Union{Label,Nothing} = nothing)
	coeftable2 = coeftbl

	# estimates
	coeftable2.cols[1] = exp.(coeftable2.cols[1])

	# standard errors
	coeftable2.cols[2] = coeftable2.cols[1] .* coeftable2.cols[2]

	# rename column1 to OR
	coeftable2.colnms[1] = "OR"

	# parse the row names and change variable names and values
	for i in 2:length(coeftable2.rownms)
		# parse row name and split into a tuple (varname, value)
		if contains(coeftable2.rownms[i],": ")
			(varname,value) = split(coeftable2.rownms[i],": ")
		else
			(varname,value) = (coeftable2.rownms[i],nothing)
		end

        # get variable label from the label dictionary
		varlabel = varlab(labels,Symbol(varname))
		vlabel = vallab(labels,Symbol(varname),parse(Int,value))

		# If value is 1 and value label is Yes, it is a binary variable
		# do not print
		if value == 1 && ismatch(r"^ *yes *$"i,vlabel)
			coeftable2.rownms[i] = varlabel
		else
			coeftable2.rownms[i] = string(varlabel, ": ", vlabel)
		end
	end

	return coeftable2
end


#--------------------------------------------------------------------------
# pairwise correlations
#--------------------------------------------------------------------------
struct pwcorr_return
    r::AbstractArray
    pval::AbstractArray
    N::AbstractArray
    colnms::Vector
end

"""
    pwcorr(a::AbstractArray)
    pwcorr(a::AbstractVector...)
    pwcorr(df::DataFrame, args::Symbol...)

Produces a n x n matrix that contain correlation coefficients between all `n` columns.
When the option `p = true` is specified, it produces a second n x n matrix containing
p-values which is calculated using the Fisher transformation. When the option `out = true`
is set, pwcorr prints the correlation coeffients and p-values (if specified) to the console
in addition to returning them in arrays.
"""
function pwcorr(a::AbstractArray)

    # if typeof(a) == DataFrame
    #     a = a[completecases(a),:]
    #     r = cor(Array(a))
    #     colnames = string.(names(a))
    # else
    #     r = cor(a)
    #     colnames = string.(1:size(a,2))
    # end
    #
    cols = size(a,2)
    N = zeros(Int64,cols,cols)
    r = zeros(Float64,cols,cols)
    pval = zeros(Float64,cols,cols)
    if typeof(a) == DataFrame
        colnames = names(a)
    else
        colnames = collect(1:size(a,2))
    end

    for j = 1:cols
        for i = j:cols
            if i != j
                if typeof(a) == DataFrame
                    x = Array(a[completecases(a[[i,j]]),[i,j]])
                else
                    x = a[:,[i,j]]
                end
                r[i,j] = r[j,i] = cor(x)[1,2]
                rows = N[i,j] = N[j,i] = size(x,1)
                z = (sqrt(rows - 3) / 2)*log((1+r[i,j])/(1-r[i,j]))
                pval[i,j] = pval[j,i] = Distributions.ccdf(Distributions.Normal(),z) / 2.0
            end
        end
    end
    return pwcorr_return(r, pval, N, colnames)
end
pwcorr(a::AbstractArray...) = pwcorr(hcat(a...))
pwcorr(a::DataFrame, args::Vector{Symbol}; out=true) = pwcorr(df[args], out = out)




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
#     df2 = df[completecases(df[[var1,var2]]),[var1,var2]]
#
#     # difference
#     d = df[var1] .- df[var2]
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
  if contains(str,"*") || contains(str,"?")
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
  maxlen = getmaxlength(vstr)

  for i=1:length(vstr)
    println(append_spaces(vstr[i],maxlen),"  ", datasize(stat(vstr[i]).size))
  end
end
