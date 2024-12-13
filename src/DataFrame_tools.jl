
#############################################################################
#
# DataFrame tools
#
#############################################################################
"""
    dfcompress(df::DataFrame)

Reduce `df`'s memory use by changing the eltype of each column in the `df` to
the type that can accommodate the larget and the smallest values within the same
integer or float class.
"""
function dfcompress(df::DataFrame; verbose = false)

    df2 = deepcopy(df)

    for v in propertynames(df)

        # if Array is empty after all missings are dropped
        # drop it from the df
        if all(ismissing.(df2[!,v]))
            if verbose
                select!(df2,Not(v))
                println(v, " was empty, now deleted")
            end                
            continue
        end

        # get the original eltype
        eltype_old = nonmissingtype(eltype(df2[!,v]))

        # if string, continue
        if eltype_old == String
            continue
        end

        # compress
        df2[!,v] = Stella.acompress(df2[!,v])

        if verbose && eltype_old != nonmissingtype(eltype(df2[!,v]))
            println(v, " was ", eltype_old, ", now ", nonmissingtype(eltype(df2[!,v])))
        end
    end

    return df2
end

"""
    acompress(da::AbstractVector)

Compresses a vector to a smallest numeric type that can hold without loss of information.
"""
function acompress(da::AbstractVector)

    # get the original eltype
    eltyp = eltype(da)
    eltype_old = nonmissingtype(eltyp)
    
    nomiss = true
    if eltyp == Union{Missing,eltype_old} && nmissing(da) > 0
        nomiss = false
    end

    # string variable - do not compress
    if eltype_old == String
        return da
    end

    if  eltype_old <: Integer

        varmin = minimum(skipmissing(da))
        varmax = maximum(skipmissing(da))
        if eltype_old != Int8 && varmin >= typemin(Int8) && varmax <= typemax(Int8)
            if nomiss
                return convert(Vector{Int8},da)
            else
                return convert(Vector{Union{Missing,Int8}},da)
            end
        end
        if eltype_old != Int16 && varmin >= typemin(Int16) && varmax <= typemax(Int16)
            if nomiss
                return convert(Vector{Int16}, da)
            else
                return convert(Vector{Union{Missing,Int16}}, da)
            end
        end
        if eltype_old != Int32 && varmin >= typemin(Int32) && varmax <= typemax(Int32)
            if nomiss
                return convert(Vector{Int32}, da)
            else
                return convert(Vector{Union{Missing,Int32}}, da)
            end
        end
    elseif eltype_old <: AbstractFloat
        # first test if the floats are actually integer numbers (all decimals are zeros)
        if all([x - floor(x) > 0 ? false : true for x in skipmissing(da)] .== true)
            if nomiss
                return acompress(convert(Vector{Int64},da))
            else
                return acompress(convert(Vector{Union{Missing,Int64}},da))
            end
        end
    end
    return da
end


function atype(df::DataFrame,v::Symbol)
    # Array type = DA for DataArray, CA for Categorical Array, and UV for Union Vector
    if isdefined(Main,:CategoricalArrays) && typeof(df[!,v]) <: CategoricalArray
        return string("CA (", replace(string(eltype(df[!,v].refs)),"UInt" => ""), ")")
    # elseif isdefined(Main,:DataArrays) && typeof(df[!,v]) <: DataArray
    #      return "DA"
    elseif isdefined(Main,:PooledArrays) && typeof(df[!,v]) <: PooledArray
         return string("PA (", replace(string(eltype(df[!,v].refs)),"UInt" => ""), ")")
    # elseif isa(eltype(df[!,v]),Union)
    #     return "UV" # Union Vector
    end
    return missing
end

function etype(df::DataFrame,v::Symbol)
    # Eltype
    if typeof(df[!,v]) <: CategoricalArray
        eltyp = string(eltype(df[!,v].pool.levels))
        if in(eltyp,["String","AbstractString"])
            eltyp = string("Str",getmaxwidth(df[!,v].pool.levels))
        end
    else
        eltyp = string(nonmissingtype(eltype(df[!,v])))
        if in(eltyp,["String","AbstractString"])
            eltyp = string("Str",getmaxwidth(df[!,v]))
        elseif eltyp == "Dates.Date"
            eltyp = "Date"
        elseif eltyp == "Dates.DateTime"
            eltyp = "DateTime"
        end
    end

    return eltyp
end

function eltype2(df::DataFrame,v::Symbol)
    if typeof(df[!,v]) <: CategoricalArray
        return eltype(df[!,v].pool.index)
    end
    return nonmissingtype(eltype(df[!,v]))
end

"""
    descr(df::DataFrame,varnames::Symbol...; nmiss::Bool = false)

Displays variables in a dataframe much like `showcols`. It can display variable labels and value labels.
It mimics Stata's `describe` command. 
"""
function descr(df::DataFrame,varnames::Symbol...; nmiss::Bool = true, max_varlen = 15, max_descr = 40)

    cutlen(str, len) = (length(str) > len ? string(str[1:len - 6],"~", str[len - 3:end]) : str)

    # get variable names
    varnames = propertynames(df)
    if length(varnames) == 0
        error("No variables in the input dataframe\n")
    end

    # output dataframe
    dfv = DataFrame(Variable = cutlen.(string.(varnames), max_varlen))
    dfv[!,:Eltype] = Vector{String}(undef,size(dfv,1))
    if nmiss
    	dfv[!,:Missing] = Vector{String}(undef,size(dfv,1))
    end

    dfv[!,:Description] = cutlen.(labels(df),max_descr)

    for (i,v) in enumerate(varnames)

        # Eltype
        dfv[i,:Eltype] = etype(df,v)

        # percent missing
        if nmiss
            _nmiss = nmissing(df[!,v])
            dfv[i,:Missing] = string(round(100 * _nmiss/size(df,1),digits=1),"%")
        end

    end

    header = ["Variable", "Eltype"]
    alignment = [:l,:l]

    if nmiss
	    header = vcat(header,"% Miss")
	    alignment = vcat(alignment,:r)
    end

    header = vcat(header,"Description")
	alignment = vcat(alignment,:l)

    # if dfout 
    # 	return dfv
    # else
        # if datalab(df) != nothing
        #     println(datalab(df))
        # end
        println("Number of observations: ", @sprintf("%12.0f",nrow(df)))
        println("Number of variables:    ", @sprintf("%12.0f",ncol(df)))
        
        pretty_table(dfv,
            alignment=alignment,
            header=header,
            crop=:none,
            vlines = [1],
            formatters = (ft_nomissing),
            show_row_number = true,
            row_number_column_title = "Column")
    # end
end

function nmissing(s::AbstractArray)
    return sum(ismissing.(s))
end

function getmaxwidth(s::AbstractArray)
    if isa(s, CategoricalArray) && nonmissingtype(eltype(s)) <: CategoricalString
	    return maximum(length.(codeunits(s.pool.levels)))
    end
	
    if nmissing(s) == size(s,1)
	    return 0
    end
	
    return  maximum(length.(codeunits.(collect(skipmissing(s)))))
end

"""
    ds(df::DataFrame, typ::Type, args...)
    ds(df::DataFrame, re::Regex)

Returns an array of Symbols of the columns whose eltype matches the `typ` or the regular expression `re`.
For `String` variables, up to two additional arguments can be specified to indicate the minimum and maximum widths
to be returned. In Julia, string variables does not have an associated widths as in Stata.
So, in this function, the width is measured by the length of the longest string. The function follows `Number` type
structure. For example, to list all numeric variables, use `Number`; for all integer variables, use `Integer`; for
all floating point numbers, use `AbstractFloat`. This function is useful to subset a DataFrame using
eltypes or column names. This function emulates the Stata `ds` command.

# Example
```
julia> showcols(aapl)
6081×7 DataFrames.DataFrame
│ Col # │ Name      │ Eltype  │ Missing │
├───────┼───────────┼─────────┼─────────┤
│ 1     │ Date      │ String  │ 0       │
│ 2     │ Open      │ Float64 │ 0       │
│ 3     │ High      │ Float64 │ 0       │
│ 4     │ Low       │ Float64 │ 0       │
│ 5     │ Close     │ Float64 │ 0       │
│ 6     │ Volume    │ Int64   │ 0       │
│ 7     │ Adj_Close │ Float64 │ 0       │

julia> ds(aapl,Int64)
1-element Array{Symbol,1}:
 :Volume

julia> ds(aapl,Float64)
5-element Array{Symbol,1}:
 :Open
 :High
 :Low
 :Close
 :Adj_Close

julia> ds(aapl,Number)
6-element Array{Symbol,1}:
 :Open
 :High
 :Low
 :Close
 :Volume
 :Adj_Close

julia> ds(aapl,String)
1-element Array{Symbol,1}:
 :Date

julia> ds(aapl,r".*lose")
2-element Array{Symbol,1}:
 :Close
 :Adj_Close

```
"""
function ds(df::DataFrame, typ::Type, args...)

    dslist = Vector{Symbol}()

    for v in propertynames(df)
        if typ == String && nonmissingtype(eltype(df[!,v])) == String
            if length(args) == 0
                push!(dslist,v)
                continue
            end

            maxlen = getmaxwidth(df[!,v])
            if length(args) == 1
                if args[1] <= maxlen
                    push!(dslist,v)
                end
            elseif length(args) == 2
                if args[1] <= maxlen <= args[2]
                    push!(dslist,v)
                end
            end
        elseif nonmissingtype(eltype(df[!,v])) == typ
            push!(dslist,v)
        elseif (typ == Integer && nonmissingtype(eltype(df[!,v])) <: Integer) || (typ == AbstractFloat && nonmissingtype(eltype(df[!,v])) <: AbstractFloat) || (typ == Number && nonmissingtype(eltype(df[!,v])) <: Number)
            push!(dslist,v)
        end
    end

    return dslist

end
function ds(df::DataFrame,re::Regex)
    dslist = Array{Symbol,1}()

    for v in propertynames(df)
        if occursin(re,string(v))
            push!(dslist,v)
        end
    end

    return dslist
end


"""
    pickone(df::DataFrame,groupvars::Vector{Symbol})

Create a vector that identifies one record in a `groupvars` in the `df` DataFrame.
This function creates a variable similar to the Stata command
`egen byte pickone = tag(groupvar)`. `groupvars` can be an array of Symbols or a single
Symbol.
"""
function pickone(df::DataFrame,groupvars::Vector{Symbol})
    df[!,:_____obs_____] = collect(1:size(df,1))
    done = zeros(Int8,size(df,1))
    for subdf in groupby(df, groupvars)
        done[subdf[1,:_____obs_____]]=1
    end
    select!(df,Not(:_____obs_____))
    return done
end
pickone(df::DataFrame,groupvar::Symbol) = pickone(df, [groupvar])


"""
    duplicates(df::DataFrame, args::Symbol...; cmd::Symbol = :report)

Reports, tags, or deletes duplicates in a DataFrame for one or more columns. Use `cmd` to
request one of the actions:

- :report - for getting a frequency table for the number of duplicates
- :tag - for identifying rows with duplicate values
- :drop - for dropping all rows with duplicate values except for the first row.
"""
function duplicates(df::DataFrame, args::Union{Symbol,String}... ; cmd::Symbol = :report)

    if in(cmd,[:report,:drop,:tag]) == false
        error("`cmd = `", cmd, "` is not supported.")
    end

    # if args are not specified, use all variables
    if length(args) == 0
        args = tuple(propertynames(df)...)
    end

    gdf = groupby(df,collect(args))
    dfx = combine(gdf,nrow => :__dups__)

    if cmd == :report
        na = freqtable(dfx, :__dups__)
        setdimnames!(na,"copies",1)
        return na
    end

    # substract 1 from :x in dfx (we are reporting 0 for unique observations)
    dfx[!,:__dups__] .-= 1
    df = leftjoin(df, dfx, on = collect(args))

    if cmd == :tag
        return df[!,:__dups__]
    end

    # if cmd == :drop
    ba = [x == 1 ? true : false for x in pickone(df,collect(args))]
    return select(df[ba,:],Not(:__dups__))
end

"""
    sample(df::AbstractDataFrame,num::Real)

Creates an Int8 array that identifies a randomly selected sample (1 = sample; 0 otherwise) from the input dataframe. `num`
specifies the sample size. If `num` is an integer, it indicates the number of rows to be selected.
If it is a float, it indicates the percentage of the input dataframe to be randomly selected.
Use `Random.seed!()` to set a seed before using `sample()`.

##Example
To select 100 rows, use

```
julia> df[!,:sample] = sample(df,100)
```

To select 20% of the original dataframe, use

```
julia> df[!,:sample2] = sample(df,.2)
```

"""
function sample(df::AbstractDataFrame,n::Real)
    len = size(df,1)
    asamp = zeros(Int8,len)
    if isa(n,AbstractFloat)
        asamp[randsubseq(collect(1:len),n)] .= 1
    else
        asamp[randperm(len)[1:n]] .= 1
    end
    return asamp
end

"""
    dfsample(df::AbstractDataFrame,num::Real)

Creates a dataframe with a randomly selected sample from the input dataframe. `num`
specifies the sample size. If `num` is an integer, it indicates the number of rows to be selected.
If it is a float, it indicates the percentage of the input dataframe to be randomly selected.
Use `Random.seed!()` to set a seed before running `dfsample()`.

##Example
To select 100 rows, use

```
julia> df2 = dfsample(df,100)
```

To select 20% of the original dataframe, use

```
julia> df2 = dfsample(df,.2)
```

"""
function dfsample(df::AbstractDataFrame,num::Real)
    return df[sample(df,num) .== 1,:]
end

"""
    dfmerge(df1::DataFrame,df2::DataFrame,linkers::Union{Symbol,Vector}; kind::Symbol = :outer)

Produces a merged dataframe with a `:_merge` variable indicating how the merge was done.
This function is a wrapper for `join` in the DataFrames package and mimics `merge ... using ...`
command in Stata. Currently, it supports a merge of two sources only.

The default merge is a `outer` join that keeps records from both sources in the merged
dataframe. `:_merge` values indicates:

- `1`: in the first source alone
- `2`: in the second source alone
- `3`: in both sources. A successful merge was achieved.

It is the responsibility of the user to keep or discard records with `:_merge` values 1 or 2.

The third argument can be either a single symbol or a vector of symbols to be used as the
linkage variables.

Optionally, `kind` argument can be specified. Available options are
`:left`, `:right`, `:inner`,`:outer`,`:semi`, and `:anti`. Consult the documentation for
the `join` function in the DataFrames package.
"""
function dfmerge(df1::DataFrame,df2::DataFrame,linkers::Union{Symbol,Vector})
    if in(:_merge,propertynames(df1))
        error("`:_merge' exists in the first dataframe")
    end
    if in(:_merge,propertynames(df2))
        error("`:_merge' exists in the second dataframe")
    end

    df1[!,:___mergeleft___] = ones(Int8,size(df1,1))
    df2[!,:___mergeright___] = ones(Int8,size(df2,1))

    df_merged = outerjoin(df1,df2,on = linkers)
    df_merged[!,:_merge] = zeros(Int8,size(df_merged,1))
    for i = 1:size(df_merged,1)
        if ismissing(df_merged[i,:___mergeright___])
            df_merged[i,:_merge] = 1
        elseif ismissing(df_merged[i,:___mergeleft___])
            df_merged[i,:_merge] = 2
        elseif df_merged[i,:___mergeleft___] == 1 && df_merged[i,:___mergeright___] == 1
            df_merged[i,:_merge] = 3
        end
    end
    select!(df_merged,Not([:___mergeleft___,:___mergeright___]))
    return df_merged
end

function get_class(val::Real,thresholds::Vector,lower::Bool)
    if lower == false
        for i = 1:length(thresholds)
            if i == 1 && val < thresholds[1]
                return 1
            elseif i > 1 && thresholds[i-1] <= val < thresholds[i]
                return i
            end
        end
    else
        for i = 1:length(thresholds)
            if i == 1 && val <= thresholds[1]
                return 1
            elseif i > 1 && thresholds[i-1] < val <= thresholds[i]
                return i
            end
        end
    end
    return length(thresholds)+1
end

"""
    classify(da::AbstractVector,thresholds::Vector; lower::Bool = false)

Produces a categorical variable based on a data array and a vector of thresholds.
Use `lower = true` to classify the threshold value to the lower class.
"""
function classify(da::AbstractVector,thresholds::Vector; lower::Bool = false)
    da2 = Vector{Union{Missing,Int8}}(undef,length(da))
    if length(thresholds) > 100
        error("Cannot use more than 100 threshold values.")
    end

    for i = 1:size(da,1)
        if ismissing(da[i])
            da2[i] = missing
        else
            da2[i] = get_class(da[i],thresholds,lower)
        end
    end
    return da2
end

"""
    addterms(fm::Formula,v::Vector{Symbol})

Adds covariates to an existing formula. `fm` is an object of `Formula` type
created by `@formula` macro. `v` is a vector of Symbols. This function is useful to
build hierarchically nested models.

### Example

```
julia> using DataFrames

julia> fm = @formula(income ~ age + male)
Formula: income ~ age + male

julia> fm2 = addterms(fm,[:race, :educ])
Formula: income ~ age + male + race + educ
```

"""
function addterms(fmm::Formula,v::Vector{Symbol})

    # create a new Formula object
    fm = deepcopy(fmm)

    # if fm.rhs is a null modle with rhs == 1 or
    # a unadjusted model with only one rhs variable,
    # convert it to an Expr first
    if fm.rhs == 1
        fm.rhs = Expr(:call,:+,Tuple(v)...)
    elseif typeof(fm.rhs) == Symbol
        fm.rhs = Expr(:call,:+,fm.rhs,Tuple(v)...)
    else
        push!(fm.rhs.args,Tuple(v)...)
    end

    return fm
end

"""
    Stella.formula(o::Symbol,v::Vector{Symbol})

Creates a Formula from a vector of covariates. This function is not exported for potential conflict with other packages.

- o - A symbol for the dependent variable
- v - A vector of symbols for covariates

## Usage

julia> fm = Stella.formula(:readmit, [:age,:male,:race,:chd,:exercise])
Formula: readmit ~ age + male + race + chd + exercise
"""
function formula(o::Symbol,v::Vector{Symbol})
    if length(v) == 1
        rhs = v[1]
    else
        rhs = Expr(:call,:+,Tuple(v)...)
    end
    e = :($o ~ $rhs)
    return Formula(e,e,o,rhs)
end

"""
    renvars!(df::DataFrame; vars::Array{Symbol,1}, case = "lower")

Rename column names in `vars` to either upper or lower cases. The default is to convert
all columns to lower case names.
"""
function renvars!(df::DataFrame; vars=[], case="lower")
    numvar = length(vars)
    symnames = propertynames(df)

    if numvar == 0
        varnames = propertynames(df)
    else
        varnames = propertynames(df[vars])
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
            rename!(df,nm => Symbol(newname))
        end
    end
end

"""
	ci(df::AbstractDataFrame,colname::Union{Symbol,String})
	ci(df::AbstractDataFrame, colname1::Union{Symbol,String}, colname2::Union{Symbol,String})
	ci(df::AbstractDataFrame,colre::Regex)

Produces column index using DataFrames.columnindex() function. It can produce a single index,
or a range of indices, or multiple indices matching a regular expression.
"""
function ci(df::AbstractDataFrame,colname::Union{Symbol,String})
    return columnindex(df,colname)
end
function ci(df::AbstractDataFrame, colname1::Union{Symbol,String}, colname2::Union{Symbol,String})
    ci1 = columnindex(df,colname1)
    ci2 = columnindex(df,colname2)
    return ci1 < ci2 ? range(ci1,ci2) : range(ci2,ci1)
end
function ci(df::AbstractDataFrame,colre::Regex)

    vec = Int[]
    for n in names(df)
        m = match(colre, n)
        if m !== nothing
            push!(vec,columnindex(df,n))
        end
    end

    return vec
end


"""
    destring(da::AbstractArray;force=true)
    destring(df::DataFrmae,strvar::Symbol;force=true)
    destring!(df::DataFrame,strvars; newvars::Vector{Symbol} = [], force=true, replace=false)

Convert a string vector to a numeric vector. Use `force = true` to coerce conversion of alphanumeric strings to
missing values. If `force = false`, any vector containing non-numeric characters are not converted.
Use `replace = true` in `destring!` to replace the original string vector with a new converted numeric vector.
If `replace` option is specified, `newvars` array is ignored.
"""
function destring(da::AbstractArray; force=true)
    if length(da) == 0
        error("Input array is empty!")
    end
    if nonmissingtype(eltype(da)) <: Number
        error("Input array is a numeric vector. No conversion needed.")
    end

    # check if the values include any alphabetic characters or decimals
    isfloat = false
    alpha = false
    da_safe = collect(skipmissing(da))
    for i in length(da_safe)
        if sum([isdigit(x) == false for x in da_safe[i]]) > 0
            alpha = true
        end
        if occursin(r"[,0-9]*\.[0-9]+",da_safe[i])
            isfloat = true
        end
    end

    if alpha && force == false
        error("Input array contains alphabetic letters. Use 'force=true' option to coerce conversion.")
    end

    T = isfloat ? Float64 : Int64
    da2 = [ismissing(x) || x == "" ? missing : parse(T,x) for x in da]

    return acompress(da2)
end
destring(df::AbstractDataFrame,strvar::Symbol; force=true) = destring(df[:,strvar],force=force)
function destring!(df::DataFrame,strvars; newvars::Vector{Symbol} = [], force=true, replace=false)

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
    rowstat(df::DataFrame,func::Function)

Creates a vector that contains the row statistic of all values on the same row of `df`
produced by the `func` function. If one of the columns contain a missing value on a row, the returned value will be mssing for that
row. This function emulates Stata's `egen` row functions such as `rowtotal`, `rowmean`, etc.

```
julia>df[:rowmean] = rowstat(df[[:var1,:var2,:var3]],mean)
```
"""
function rowstat(df::DataFrame,func::Function)

    # df size
    nrow, ncol = size(df)

    # array to be returned
    da = zeros(Union{Missing,Float64},size(df,1))

    # temporary array of row values
    ta = Vector{Union{Missing,Float64}}(undef,size(df,2))

    # iterate over all rows
    for i = 1:nrow

        # counter for non-missing values in each row
        k=0
        for j = 1:ncol
            if ismissing(df[i,j]) == false
                k += 1
                ta[k] = df[i,j]
            end
        end

        # if any of the columns has a missing value, return missing
        if k < ncol
            da[i] = missing
        else
            tmpfloat = func(ta[1:k])
            da[i] = isnan(tmpfloat) ? missing : tmpfloat
        end
    end
    return da
end

"""
    xtile(da::AbstractArray;nq::Int = 4, lower::Bool = false, cutoffs::Union{Nothing,AbstractVector} = nothing)
    xtile(df::DataFrame,varname::Symbol;nq::Int = 4, lower::Bool = false, cutoffs::Union{Nothing,AbstractVector} = nothing)

Create a DataArray{Int8,1} that identifies `nq` categories based on values in `da`.
The default `nq` is 4. `cutoffs` vector can be provided to make custom categories.
`cutoffs` vector is expected to contain `nq - 1` elements. The minimum and maximum values
will be computed. `lower = true` will make the threshold values to go in the lower category (default: false).

```
julia> df[:agecat] = xtile(df[:age], nq = 3)
```
"""
function xtile(da::AbstractArray ; nq::Int = 4, lower::Bool = false, cutoffs::Union{Nothing,AbstractVector} = nothing)

    if nq == 1
        error("`nq` must be greater than 1")
    end

    da2 = collect(skipmissing(da))
    if cutoffs == nothing
	    cutoffs = nquantile(da2,nq)
    elseif length(cutoffs) == nq - 1
        cutoffs = vcat(minimum(da2), cutoffs, maximum(da2))
    else
        error("`cutoffs` vector length is not consistent with `nq`. It must be 1 greater or 1 less than `nq`.")
    end

    if lower
        return [ismissing(x) ? missing : qval_low(x,cutoffs) for x in da2]
    else
        return [ismissing(x) ? missing : qval_high(x,cutoffs) for x in da2]
    end
end
xtile(df::DataFrame,arg::Symbol; nq::Int = 4, cutoffs::Union{Nothing,AbstractVector} = nothing) = xtile(df[!,arg], nq = nq, cutoffs = cutoffs)

function qval_high(val::Real,cut::Vector)
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

function qval_low(val::Real,cut::Vector)
    cl = length(cut)

    for i in 2:cl-1
        if i == 2 && val <= cut[i]
            return 1
        elseif i == cl-1 && cut[i] <= val
            return i
        elseif cut[i] < val <= cut[i+1]
            return i
        end
    end
    warn("Error - check qval function")
end

"""
	categorical!(df::AbstractDataFrame,vv::Union{Symbol,Vector{Symbol}})

converts vectors into CategoricalArrays
"""
function categorical!(df::AbstractDataFrame,vv::Union{Symbol,Vector{Symbol}})
    for v in vcat(vv)
        if isa(df[:,v],CategoricalArray) == false
            df[!,v] = categorical(df[!,v])
        end
    end
end


"""
	uncategorical!(df::AbstractDataFrame,vv::Union{Symbol,Vector{Symbol}})

uncategorizes a vector of CategoricalArrays by replacing the refs with their original values.
"""
function uncategorical!(df::AbstractDataFrame,vv::Union{Symbol,Vector{Symbol}})
    for v in vcat(vv)
        if isa(df[:,v],CategoricalArray)
            lev = df[:,v].pool.levels
            df[!,v] = [x == 0 ? missing : lev[x] for x in df[:,v].refs]
        end
    end
end

"""
	uncategorize(v::AbstractArray)

produces an array by uncategorizing a CategoricalArray by replacing the refs with its original values.
"""
function uncategorize(v::CategoricalArray)
    return [x == 0 ? missing : v.pool.levels[x] for x in v.refs]
end

"""
    rowfirst(df::AbstractDataFrame,vars::AbstractArray)
    rowlast(df::AbstractDataFrame,vars::AbstractArray)
    rowmax(df::AbstractDataFrame,vars::AbstractArray)
    rowmin(df::AbstractDataFrame,vars::AbstractArray)
    rowmean(df::AbstractDataFrame,vars::AbstractArray)
    rowmedian(df::AbstractDataFrame,vars::AbstractArray)
    rowmiss(df::AbstractDataFrame,vars::AbstractArray)
    rownonmiss(df::AbstractDataFrame,vars::AbstractArray)
    rowpctile(df::AbstractDataFrame,vars::AbstractArray,pct )
    rowsd(df::AbstractDataFrame,vars::AbstractArray)
	rowtotal(df::AbstractDataFrame,vars::AbstractArray)

produces an array that contains the statistics based on all variables in the AbstractArray, 
ignoring missing values. This is a Julia implementation of Stata "row" egen functions.
"""
function rowfirst(df::AbstractDataFrame,vars::AbstractArray)
    return [collect(skipmissing(x))[1] for x in eachrow(df[:,vars])]
end

function rowlast(df::AbstractDataFrame, vars::AbstractArray)
    return [collect(skipmissing(x))[end] for x in eachrow(df[:, vars])]
end

function rowmax(df::AbstractDataFrame, vars::AbstractArray)
    return [maximum(collect(skipmissing(x))) for x in eachrow(df[:, vars])]
end

function rowmin(df::AbstractDataFrame, vars::AbstractArray)
    return [minimum(collect(skipmissing(x))) for x in eachrow(df[:, vars])]
end

function rowmean(df::AbstractDataFrame, vars::AbstractArray)
    return [mean(collect(skipmissing(x))) for x in eachrow(df[:, vars])]
end

function rowmedian(df::AbstractDataFrame, vars::AbstractArray)
    return [median(collect(skipmissing(x))) for x in eachrow(df[:, vars])]
end

function rowmiss(df::AbstractDataFrame, vars::AbstractArray)
    return sum.(ismissing.(eachrow(df[:, vars])))
end

function rownonmiss(df::AbstractDataFrame, vars::AbstractArray)
    return length(vars) .- sum.(ismissing.(eachrow(df[:, vars])))
end

function rowsd(df::AbstractDataFrame, vars::AbstractArray)
    return [std(collect(skipmissing(x))) for x in eachrow(df[:, vars])]
end

function rowpctile(df::AbstractDataFrame, vars::AbstractArray, p::Float64 = 0.5)
    return [quantile(collect(skipmissing(x)), p) for x in eachrow(df[:, vars])]
end

function rowtotal(df::AbstractDataFrame, vars::AbstractArray)
    return [sum(collect(skipmissing(x))) for x in eachrow(df[:, vars])]
end

function firstrow(df::AbstractDataFrame, groupid::Symbol)
    keep = falses(nrow(df))
    for i = 1:nrow(df)
        if i == 1 || df[i, groupid] != df[i-1, groupid]
            keep[i] = true
        end
    end
    return df[keep.==true, :]
end
function firstrow(df::AbstractDataFrame, groupids::Vector{Symbol})
    keep = falses(nrow(df))
    diff = falses(length(groupids))
    for i = 1:nrow(df)
        if i == 1
            keep[i] = true
            continue
        end

        diff .= false
        for gid in groupids
            if df[i, gid] != df[i-1, gid]
                diff[j] = true
            end
        end

        if all(diff)
            keep[i] = true
        end
    end

    return df[keep.==true, :]
end


