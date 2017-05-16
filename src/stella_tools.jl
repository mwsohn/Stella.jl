using DataFrames, FreqTables

"""
    prepend_spaces(str,maxlength)

Create a string of length `maxlength`. When `str` is shorter than `maxlength`,
blank spaces are prepended. When `str` is longer than `maxlength`, it is
truncated to fit within the `maxlength`. This function is used to right-justify
`str` in output.

# Example
```jldoctest
julia> prepend_spaces("test",8)
"    test"

julia> prepend_spaces("test result",8)
"test res"
```
"""
function prepend_spaces(str,maxlength)
    len = length(str)
    if maxlength < len
        return str[1:maxlength]
    end
    return string(repeat(" ",maxlength-len),str)
end


"""
    append_spaces(str,maxlength)

Create a string of length `maxlength`. When `str` is shorter than `maxlength`,
blank spaces are appended at the end. When `str` is longer than `maxlength`, it is
truncated to fit within the `maxlength`. This function is used to left-justify
`str` in output.

# Example
```jldoctest
julia> append_spaces("test",8)
"test    "

julia> append_spaces("test result",8)
"test res"
```
"""
function append_spaces(str,maxlength)
    len = length(str)
    if maxlength < len
        return str[1:maxlength]
    end
    return string(str,repeat(" ",maxlength-len))
end


"""
    smallest(da::DataArray, n::Int) or smallest(df::DataFrame, varname::Symbol, n::Int)

List the `n` smallest values in `da` in a descending order. The default `n` is 5.
"""
function smallest(da::DataArray; n::Int = 5)
	return sort(dropna(da))[1:n]
end
smallest(df::DataFrame,varname::Symbol; n = 5) = smallest(df[varname], n = n)

"""
    largest(da::DataArray, n::Int) or largest(df::DataFrame, varname::Symbol, n::Int)

List the `n` largest values in `da` in an ascending order. The default `n` is 5.
"""
function largest(da::DataArray; n::Int = 5)
	return sort(dropna(da))[end-n+1:end]
end
largest(df::DataFrame,varname::Symbol; n = 5) = largest(df[varname], n = n)

"""
    univariate(da::DataArray) or univariate(df::DataFrame,varname::Symbol)

Produce univariate statistics from the `da` and return a DataFrame with two columns,
:Statistic and :Value. Computed statistics include `N_Total` (number of rows
in the input DataArray), `N_Missing` (number of rows with NA's), `N_Used` (number of non-NA rows),
`Sum` (sum), `Mean` (mean), `SD` (standard deviation), `Var` (variance), `Min` (minimum),
`P25` (25th percentile), `Median` (median), `P75` (75th percentile), `Max` (maximum),
`Skewness` (skewness), and `Kurtosis` (kurtosis).
"""
function univariate(da::AbstractArray)
	n_all = size(da,1)
	da2 = da[!isna(da)]
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
		Value = [n_all,
			n_all - size(da2,1),
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



"""
    tab(x::AbstractArray...; rmna::Bool = true, weights::AbstractVector = UnitWeights())
    tab(df::DataFrame,vars::Symbol...; rmna::Bool = true, weights::AbstractVector = UnitWeights())

Produce n-way frequency table from a DataFrame or any type of arrays. Weights can be used to obtain
weighted frequencies. `tab` is mainly a wrapper for the excellent `FreqTables` package for all
but DataArrays with NA values. Use `rmna = false` to obtain frequencies that include NA values.
The returned table is in `NamedArrays`. Frequencies are in an n-dimensional array `na.array`
where `na` is the returned NamedArray. `na.dimnames` contain row and column values.
When NA values are included, the `dimnames` (see `NamedArrays` package) are returned in string
values rather than in the original data type.
"""
function tab(df::DataFrame,args::Symbol...;rmna = true, weights::AbstractVector = FreqTables.UnitWeights())
    if rmna == true
        df = df[completecases(df[collect(args)]),collect(args)]
    end

    # find out if the args columns contain any NA values
    cc = rmna ? trues(size(df,1)) : completecases(df[collect(args)])
    if rmna || sum(cc) == size(df,1)
        a = freqtable([df[y] for y in args]...; weights = weights)
        setdimnames!(a, args)
        return a
    end

    # there are NA values and so we cannot use freqtable
    a = tab([df[y] for y in args]...) #; weights = weights)
    setdimnames!(a,args)
    return a
end


function tabv2(df::DataFrame,args::Symbol...;rmna = true, weights::AbstractVector = FreqTables.UnitWeights())
    if rmna == true
        df = df[completecases(df[collect(args)]),collect(args)]
    end

    # find out if the args columns contain any NA values
    cc = rmna ? trues(size(df,1)) : completecases(df[collect(args)])
    if rmna || sum(cc) == size(df,1)
        a = freqtable([df[y] for y in args]...; weights = weights)
        setdimnames!(a, args)
        return a
    end

    # there are NA values and so we cannot use freqtable
    # split the data into two and combine them later
    df1 = df[cc,:]
    df2 = df[cc .== false,:]
    a = freqtable([df[y] for y in args]...)
    b = __tab([df[y] for y in args]...)

    # combine them
    # indices
    a_dimnames = [string.(x) for x in names(a)]



    setdimnames!(a,args)
    return a
end

tab(args::AbstractVector...; weights::AbstractVector = FreqTables.UnitWeights()) = ___tab(args) #;weights=weights)
function ___tab(x::NTuple) #; weights::AbstractVector = FreqTables.UnitWeights())
    # if isa(weights,UnitWeights) == false
    #     error("Frequencies with NA values cannot be weighted.")
    # end

    ncols = length(x)
    nrows = length(x[1])
    for i = 2:ncols
        if nrows != length(x[i])
            error("Columns do not have the same rows.")
        end
    end

    # create output arrays
    vdims = Array{Array,1}(ncols)
    vnums = zeros(Int,ncols)
    for i = 1:ncols
        if sum(x[i].na) > 0 # there are NA values in this vector
            vdims[i] = sort(levels(dropna(x[i])))
            vnums[i] = size(vdims[i],1) + 1 # one for the NA
        else
            vdims[i] = sort(levels(x[i]))
            vnums[i] = size(vdims[i],1)
        end
    end

    # allocate memory for the output array
    a = zeros(Int,vnums...)
    idxvec = zeros(Int,ncols)

    # get frequencies with the NA values in the vnums... cell for each column
    for i = 1:nrows
        for j = 1:ncols
            @inbounds idxvec[j] = isna(x[j][i]) ? vnums[j] : findfirst(vdims[j],x[j][i])
        end
        @inbounds a[idxvec...] += 1
    end

    dimnames = Vector{Array}(ncols)
    for i = 1:ncols
        dimnames[i] = Vector{String}(vnums[i])
        for j = 1:length(vdims[i])
            dimnames[i][j] = string(vdims[i][j])
        end
        if vnums[i] > length(vdims[i])
            dimnames[i][vnums[i]] = "NA"
        end
    end

    return NamedArray(a, tuple(dimnames...), ntuple(i -> "Dim$i", ncols))

end


"""
    tabstat(df::DataFrame,varname::Symbol,groupvar::Symbol)

Produce a DataFrame that contains summary statistics for each `groupvar` subgroup
of the `varname` column in the `df`. The following are computed: `n` (total non-missing rows),
`mean` (mean), `sd` (standard deviation), `min` (minimum), `p25` (25th percentile),
`median` (median), `p75` (75th percentile), and `max` (maximum).
"""
function tabstat(indf::DataFrame, var1::Symbol, groupvar::Symbol; s::Vector{Function} = [N,mean,sd,minimum,p25,median,p75,maximum ])

    if length(s) == 0
        error("No statistic function was specified.")
    end

    namevec = [Symbol(replace(string(x),"Stella.","")) for x in s]


    gvnum = length(levels(indf[groupvar]))

    outdf = DataFrame()
    outdf[groupvar] = DataArray(eltype(groupvar),gvnum)
    for j = 1:length(namevec)
        if namevec[j] == :N
            outdf[namevec[j]] = DataArray(Int,gvnum)
        else
            outdf[namevec[j]] = DataArray(Float64,gvnum)
        end
    end

    i = 1
    for subdf in groupby(indf, groupvar)
        da = dropna(subdf[var1])
        if length(da) == 0
            continue
        end
        outdf[i,groupvar] = subdf[1,groupvar]
        for j = 1:length(namevec)
            outdf[i,namevec[j]] = s[j](da)
        end
        i += 1
    end
    return outdf
end


"""
    toNA(da::DataArray,v::Array{Symbol,1})

Convert values in `v` array to NAs in the `da` DataArray and return a new DataArray.
"""
function toNA(da::DataArray, varray = [])
    if isempty(varray)
        return da
    end

    for i = 1:length(da)
        if isna(da[i])
            continue
        elseif in(da[i],varray)
            da[i] = NA
        end
    end
    return da
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
    return DataArray(done)
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

```jldoctest
julia> df[:testmean] = substat(df,:test, :class, mean)
```
"""
function substat(df::DataFrame, varname::Symbol, groupvars::Vector{Symbol}, func::Function)
    if (eltype(df[varname]) <: Number) == false
        error("Only numeric variables are allowed.")
    end

    df2 = by(df,groupvars) do subdf
        df3 = dropna(subdf[varname])
        if size(df3,1) == 0
            DataFrame(x1 = NaN)
        else
            DataFrame(x1 = func(df3))
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


#----------------------------------------------------------------------------
# stats by subdataframe - end
#----------------------------------------------------------------------------


"""
    renvars!(df::DataFrame; vars::Array{Symbol,1}, case = "lower")

Rename column names in `vars` to either upper or lower cases. The default is to convert
all columns to lower case names.
"""
function renvars!(df::DataFrame; vars=[], case="lower")
    numvar = length(vars)
    symnames = names(df)

    if numvar == 0
        varnames = names(df)
    else
        varnames = names(df[vars])
    end

    for nm in varnames

        if case == "lower" || case == "LOWER"
            newname = lowercase(string(nm))
        elseif case == "upper" || case == "UPPER"
            newname = uppercase(string(nm))
        else
            error(option," ", case," is not allowed.")
        end
        if nm != Symbol(newname)
            rename!(df,nm,Symbol(newname))
        end
    end
end

function vidx(df::DataFrame,varname::Symbol)
    varnames = names(df)
    for i=1:length(varnames)
        if varnames[i] == varname
            return i
        end
    end
    return 0 # not found
end

"""
    destring(da::DataArray;force=true)
    destring(df::DataFrmae,strvar::Symbol;force=true)
    destring!(df::DataFrame,strvars::Array{Symbol,1}; newvars::Array{Symbol,1} = [], force=true, replace=true)

Convert a string DataArray to a numeric DataArray. Use `force = true` to coerce conversion of alphanumeric strings to
`NA` values. If `force = false`, any DataArray containing non-numeric values are not converted.
Use `replace = true` in `destring!` to replace the original string DataArray with a new converted numeric DataArray.
If `replace` option is specified, `newvars` array is ignored.
"""
function destring(da::DataArray; force=true)
    if eltype(da) <: Number
        error(da," is a numeric DataArray.")
    end

    # check if the values include any alphabetic characters or decimals
    isfloat = false
    alpha = false
    for i in length(da)
        if isalpha(da[i])
            alpha = true
        end
        if ismatch(r"[,0-9]*\.[0-9]+",da[i])
            isfloat = true
        end
    end

    if alpha && force == false
        error(arg," contains alphabetic letters. Use 'force=true' option to coerce conversion.")
    end

    T = isfloat ? Float64 : Int
    da2 = DataArray(T,size(da,1))

    for i in 1:size(da,1)
        da2[i] = isna(da[i]) || isalpha(da[i]) ? NA : parse(T,da[i])
    end

    return dacompress(da2)
end
destring(df::DataFrame,strvar::Symbol; force=true) = destring(df[strvar],force=force)

function destring!(df::DataFrame,strvars::Array{Symbol,1}; newvars::Array{Symbol,1} = [], force=true, replace=false)

    if replace
        for v in strvars
            df[v] = destring(df[v],force=force)
        end
    else
        # check if there are same number of symbols in strvars
        if length(strvars) != length(newvars)
            error("The number of symbols in ", strvars, " and ", newvars, " are not the same.")
        end
        for i in 1:length(strvars)
            df[newvars[i]] = destring(df[strvars[i]], force=force)
        end
    end
end


"""
    rowsum(df::DataFrame)

Creates a DataArray that contains the row total of all values on the same row of `df`.
If one of the columns contain an NA value on a row, an NA value will be returned for that
row. This function emulates Stata's `egen rowsum = rowtotal(var1 - var3)`.

```jldoctest
julia>df[:rowsum] = df[[:var1,:var2,:var3]]
```

If the position numbers for `:var1` (e.g., 4), `:var2` (5), `:var3` (6) are known and consecutive,
you can specify them as follows:

```jldoctest
julia>df[:rowsum] = df[collect(4:6)]
```
"""
function rowsum(df::DataFrame)

    isfloat = false
    for i in 1:size(df,2)
        if eltype(df[i]) <: AbstractFloat
            isfloat = true
        end
    end

    if isfloat
        da = DataArray(zeros(Float64,size(df,1)))
    else
        da = DataArray(zeros(Int64,size(df,1)))
    end

    ba = complete_cases(df)
    for i = 1:size(df,1)
        if ba[i] == false
            da[i] = NA
        else
            for j = 1:size(df,2)
                da[i] += df[i,j]
            end
        end
    end
    return da
end

"""
    rowstat(df::DataFrame,func::Function)

Creates a DataArray that contains the row statistic of all values on the same row of `df`
produced by the `func` function. If one of the columns contain an NA value on a row, an NA value will be returned for that
row. This function emulates Stata's `egen` row functions such as `rowtotal`, `rowmean`, etc.

```jldoctest
julia>df[:rowmean] = rowstat(df[[:var1,:var2,:var3]],mean)
```

If the position numbers for `:var1` (e.g., 4), `:var2` (5), `:var3` (6) are known and consecutive,
you can specify them as follows:

```jldoctest
julia>df[:rowstd] = rowstat(df[collect(4:6)],std)
```
"""
function rowstat(df::DataFrame,func::Function)

    isfloat = false
    for i in 1:size(df,2)
        if eltype(df[i]) <: AbstractFloat
            isfloat = true
        end
    end

    if isfloat
        da = DataArray(zeros(Float64,size(df,1)))
        ta = Array{Float64,1}(size(df,2))
    else
        da = DataArray(zeros(Int64,size(df,1)))
        ta = Array{Int64,1}(size(df,2))
    end

    ba = complete_cases(df)
    for i = 1:size(df,1)
        if ba[i] == false
            da[i] = NA
        else
            for j = 1:size(df,2)
                ta[j] = df[i,j]
            end
            da[i] = func(ta)
        end
    end
    return da
end

"""
    xtile(da::DataArray;nq::Int = 4, cutoffs::Union{Void,AbstractVector} = nothing)
    xtile(df::DataFrame,varname::Symbol;nq::Int = 4, cutoffs::Union{Void,AbstractVector} = nothing)

Create a DataArray{Int8,1} that identifies `nq` categories based on values in `da`.
The default `nq` is 4. `cutoffs` vector can be provided to make custom categories.
`cutoffs` vector is expected to contain `nq - 1` elements. The minimum and maximum values
will be computed.

```jldoctest
julia> df[:agecat] = xtile(df[:age], nq = 3)
```
"""
function xtile(da::DataArray ; nq::Int = 4, cutoffs::Union{Void,AbstractVector} = nothing)

	function qval(val::Real,cut::Vector)
        cl = length(cut)
		for i in 2:cl-1
            if i == 2 && val < cut[i]
                return 1
            elseif i == cl-1 && cut[i] <= val
                return i
            elseif cut[i] <= val < cut[i+1]
				return i
			end
		end
		warn("Error - check qval function")
	end
    if cutoffs == nothing
	    cutoffs = nquantile(dropna(da),nq)
    elseif length(cutoffs) == nq - 1
        cutoffs = vcat(minimum(dropna(da)), cutoffs, maximum(dropna(da)))
    else
        error("`cutoffs` vector length is not consistent with `nq`. It must be 1 greater or 1 less than `nq`.")
    end

	return convert(DataArray{Int8,1},DataArray([isna(x) ? NA : qval(x,cutoffs) for x in da]))
end
xtile(df::DataFrame,arg::Symbol; nq::Int = 4, cutoffs::Union{Void,AbstractVector} = nothing) = xtile(df[arg], nq = nq, cutoffs = cutoffs)

#
# function xtile(da::DataArray ; nq::Int = 4)
#
# 	function qval(val::Real,cut::Vector)
# 		for i in 1:length(cut)-1
# 			if cut[i] <= val <= cut[i+1]
# 				return i
# 			end
# 		end
# 		warn("Error - check qval function")
# 	end
# 	cutoffs = nquantile(dropna(da),nq)
#
# 	return convert(DataArray{Int8,1},DataArray([isna(x) ? NA : qval(x,cutoffs) for x in da]))
# end
# xtile(df::DataFrame,arg::Symbol; nq::Int = 4) = xtile(df[arg], nq = nq)
#




# recode a categorical variable to new values
function recode(da::DataArray, coding::Dict)
    #df[new] = DataArray([isna(x) ? NA : coding[x] for x in df[old]])
    # first create another dict whose keys do not include arrays
    coding2 = Dict()
    for key in keys(coding)
        for k in key
            coding2[k] = coding[key]
        end
    end

    # create a new column
    # df[new] = Array{Int8,1}(size(df,1))
    ra = Array{Any,1}(size(df,1))
    for i in 1:size(df,1)
        ra[i] = haskey(coding2,da[i]) ? coding2[da[i]] : da[i]
    end

    return convert(DataArray{Int8,1},ra)
end
recode(df::DataFrame,varname::Symbol,coding::Dict) = recode(df[varname],coding)

#----------------------------------------------------------------------------
# eform
#----------------------------------------------------------------------------
function eform(coeftbl::StatsBase.CoefTable)
	coeftable2 = coeftbl

	# estimates
	coeftable2.cols[1] = exp(coeftable2.cols[1])

	# standard errors
	coeftable2.cols[2] = coeftable2.cols[1] .* coeftable2.cols[2]

	# rename column1 to OR
	coeftable2.colnms[1] = "OR"

	return coeftable2
end


function _getval(dt::Dict,val)
    return haskey(dt,val) ? dt[val] : val
end

function eform(coeftbl::StatsBase.CoefTable, label_dict::Dict)
	coeftable2 = coeftbl

	# estimates
	coeftable2.cols[1] = exp(coeftable2.cols[1])

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
		varlabel = _getval(label_dict["variable"],varname)
        if varlabel == ""
            varlabel = varname
        end

        # get value labels
		if value == nothing
			coeftable2.rownms[i] = varlabel
			continue
		end

		lblname = haskey(label_dict["label"], varname) ? label_dict["label"][varname] : nothing

		value2 = parse(Int,value)
		vlabel = (lblname != nothing && haskey(label_dict["value"],lblname) && haskey(label_dict["value"][lblname],value2)) ?
			_getval(label_dict["value"][lblname],value2) : value2

		# If value is 1 and value label is Yes, it is a binary variable
		# do not print " - 1"
		if value2 == 1 && ismatch(r"^ *yes *$"i,vlabel)
			coeftable2.rownms[i] = varlabel
		else
			coeftable2.rownms[i] = string(varlabel, ": ", vlabel)
		end
	end

	return coeftable2
end


#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# Sum of squares
#---------------------------------------------------------------------------
function ss(da)
	tmean = mean(da)
	return mapreduce(x->(x-tmean)^2,+,da)
end

#--------------------------------------------------------------------------
# anovap - return p value only after an one-way ANOVA
#--------------------------------------------------------------------------
function anovap(df::DataFrame,dep::Symbol,cat::Symbol)
	sstotal = ss(df[dep])
	ssbetween_dataframe = by(df,cat, df2 -> ss(df2[dep]))
	ssbetween = sum(ssbetween_dataframe[:x1])
	sswithin = sstotal - ssbetween
	dfwithin = size(ssbetween_dataframe,1) - 1
	dftotal = size(df[dep],1) - 1
	dfbetween = dftotal - dfwithin
	mswithin = sswithin / dfwithin
	msbetween = ssbetween / dfbetween
	fstat = mswithin / msbetween

	return ccdf(FDist(dfwithin,dfbetween),fstat)
end

#--------------------------------------------------------------------------
# anova - return an one-way ANOVA table in a DataFrame
#--------------------------------------------------------------------------
function anova(df::DataFrame,dep::Symbol,cat::Symbol)
	sstotal = ss(df[dep])
	ssbetween_dataframe = by(df,cat, df2 -> ss(df2[dep]))
	ssbetween = sum(ssbetween_dataframe[:x1])
	sswithin = sstotal - ssbetween
	dfwithin = size(ssbetween_dataframe,1) - 1
	dftotal = size(df[dep],1) - 1
	dfbetween = dftotal - dfwithin
	mswithin = sswithin / dfwithin
	msbetween = ssbetween / dfbetween
	fstat = mswithin / msbetween

	df2 = DataFrame(
		Source = ["Between","Within","Total"],
		SS = [sswithin,ssbetween,sstotal],
		df = [dfwithin,dfbetween,dftotal],
		MS = [mswithin,msbetween,0.],
		F = [fstat,0.,0.],
		P = [ccdf(FDist(dfwithin,dfbetween),fstat),0.,0.])

	# assign NA's to empty cells
	df2[3,:MS] = NA
	df2[2:3,:F] = NA
	df2[2:3,:P] = NA

	return df2
end



#--------------------------------------------------------------------------
# pairwise correlations
#--------------------------------------------------------------------------
"""
    pwcorr(a::AbstractArray; p = false)
    pwcorr(a::AbstractArray...; p = false)
    pwcorr(df::DataFrame; p = false, out=true)
    pwcorr(df::DataFrame, args::Symbol...; p = false, out=true)

Produces a n x n matrix that contain correlation coefficients between all `n` columns.
When the option `p = true` is specified, it produces a second n x n matrix containing
p-values which is calculated using the Fisher transformation. When the option `out = true`
is set, pwcorr prints the correlation coeffients and p-values (if specified) to the console
in addition to returning them in arrays.
"""
function pwcorr(a::AbstractArray; p::Bool = false)

    r = cor(a)

    if p == false
        return r
    end

    (rows,cols) = size(a)
    pval = zeros(Float64,cols,cols)
    for j = 1:cols
        for i = j:cols
            if i != j
                z = (sqrt(rows - 3) / 2)*log((1+r[i,j])/(1-r[i,j]))
                pval[i,j] = pval[j,i] = Distributions.ccdf(Distributions.Normal(),z) / 2.0
            end
        end
    end
    return r, pval
end
pwcorr(a::AbstractArray...; p::Bool = false) = pwcorr(hcat(a...), p = p)
function pwcorr(a::DataFrame; p::Bool = false, out = true)
    b = pwcorr(Array(a[completecases[a],:]), p = p)

    if out == false
        return b
    end

    pwcorr_print(b, p, string.(names(a)))
    return b
end
pwcorr(a::DataFrame, args::Symbol...; p::Bool = false, out=true) = pwcorr(df[collect(args)], p = p, out = out)

function pwcorr_print(b::Array, p::Bool, column_names::AbstractVector{String}; width::Int = 9)

    if p == true
        bc = b[1]
        bp = b[2]
    else
        bc = b
    end

    ncol = size(bc,1)

    print(repeat(" ",width+1))

    for i = 1:ncol
        print(prepend_spaces(column_names[i],width))
        if i < ncol
            print(" ")
        end
    end
    print("\n")

    for i = 1:ncol
        print(prepend_spaces(column_names[i],width)," ")
        for j=1:i
            print(@sprintf("%9.4f",bc[i,j])," ")
            if i == j
                print("\n")
            end
        end
        if p == true
            print(repeat(" ",width+1))
            for j=1:i
                print(@sprintf("%9.4f",bp[i,j])," ")
                if i == j
                    print("\n")
                end
            end
        end
    end
end