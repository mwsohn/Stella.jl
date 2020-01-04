function memory_saved(dt,n)

    # total size and levels
    tsize = sum([ismissing(x) ? 0 : length(x.value) for x in dt.data[n]])
    s = Set(dt.data[n][:].values)

    # for refs size
    refsize = length(s) < typemax(UInt8) ? 1 : length(s) < typemax(UInt16) ? 2 : 4

    # approximate size in categorical arrays
    rsize = sum(length.(collect(s))) + refsize*dt.rows

    return tsize*.8 > rsize ? true : false
end

function readstat2dataframe(dt::ReadStat.ReadStatDataFrame,verbose=false)
    df = DataFrame()

    dateoffset = Date(1960,1,1)
    datetimeoffset = DateTime(1960,1,1,0,0,0)

    for i = 1:dt.columns
        if verbose
            println("Processing ",i," ",dt.headers[i])
        end

        # Date or DateTime values
        if dt.types[i] <: Integer && dt.formats[i] in ("%d","%td")
            df[!,dt.headers[i]] = [x.hasvalue ? dateoffset + Dates.Day(x.value) : missing for x in dt.data[i]]
        elseif dt.types[i] <: Integer && dt.formats[i] in ("%tc","%tC")
            df[!,dt.headers[i]] = [x.hasvalue ? datetimeoffset + Dates.Millisecond(x.value) : missing for x in dt.data[i]]
        elseif dt.types[i] <: AbstractString
            if memory_saved(dt,i)
                df[!,dt.headers[i]] = categorical(String[x.hasvalue ? x.value : "" for x in dt.data[i]],true)
            else
                df[!,dt.headers[i]] = [x.hasvalue || x.value != "" ? x.value : missing for x in dt.data[i]]
            end
        # elseif dt.types[i] == Int8 && dt.val_label_keys[i] != ""
        #     # variable is Int8 type and has label name, convert it to categorical array
        #     df[dt.headers[i]] = categorical(Union{Missing,dt.types[i]}[x.hasvalue ? x.value : missing for x in dt.data[i]],true)
        else
            df[!,dt.headers[i]] = Union{Missing,dt.types[i]}[x.hasvalue ? x.value : missing for x in dt.data[i]]
        end
    end
    return df
end
function rs2df(dt::ReadStat.ReadStatDataFrame)
    df=DataFrame()
    for i=1:dt.columns
        df[!,dt.headers[i]] = Union{Missing,dt.types[i]}[x.hasvalue ? x.value : missing for x in dt.data[i]]
    end
end

import DataFrames.DataFrame
DataFrame(x::ReadStat.ReadStatDataFrame) = rs2df(x)

function get_labels(dt::ReadStat.ReadStatDataFrame)
    # labels
    label = Label()

    for i=1:dt.columns
        label.var[dt.headers[i]] = dt.labels[i]
        if dt.val_label_keys[i] != ""
            label.lblname[dt.headers[i]] = Symbol(dt.val_label_keys[i])
            label.val[Symbol(dt.val_label_keys[i])] = dt.val_label_dict[dt.val_label_keys[i]]
        end
    end

    return label
end

"""
    df, label = read_stata(fn::String)

Converts a Stata datafile into a Julia DataFrame. It produces two memory objects (`df` and `label`).
The first is the dataframe. The second is a [Lables](https://github.com/mwsohn/Labels.jl) object
that provide functionality to attach variable and value labels.
"""
function read_stata(fn::String)

    dt = ReadStat.read_dta(fn)
    df = readstat2dataframe(dt)

    return df,get_labels(dt)
end


#############################################################################
#
# DataFrame tools
#
#############################################################################
"""
    dfcompress!(df::DataFrame)

Reduce `df`'s memory use by changing the eltype of each column in the `df` to
the type that can accommodate the larget and the smallest values within the same
integer or float class.
"""
function dfcompress!(df::DataFrame)
    for v in names(df)

        # if Array is empty after all missings are dropped
        # drop it from the df
        if length(df[!,v]) == sum(ismissing.(df[!,v]))
            delete!(df,v)
            println(v, " was empty, now deleted")
            continue
        end

        # get the original eltype
        eltype_old = nonmissingtype(eltype(df[!,v]))

        # if string, continue
        if eltype_old == String
            continue
        end

        # compress
        df[!,v] = acompress(df[:,v])

        if eltype_old != nonmissingtype(eltype(df[!,v]))
            println(v, " was ", eltype_old, ", now ", nonmissingtype(eltype(df[!,v])))
        end
    end
end

"""
    acompress(da::AbstractVector)

Compresses a vector to a smallest numeric type that can hold without loss of information.
"""
function acompress(da::AbstractVector)

    # get the original eltype
    eltype_old = nonmissingtype(eltype(da))

    # string variable - do not compress
    if eltype_old == String
        return da
    end

    # copy the vector
    da2 = copy(da)

    if  eltype_old <: Integer
        # get minimum and maximum values
        varmin = minimum(collect(skipmissing(da2)))
        varmax = maximum(collect(skipmissing(da2)))
        # test if Int8
        if eltype_old != Int8 && varmin >= typemin(Int8) && varmax <= typemax(Int8)
            return convert(Vector{Union{Missing,Int8}},da2)
        elseif eltype_old != Int16 && varmin >= typemin(Int16) && varmax <= typemax(Int16)
            return convert(Vector{Union{Missing,Int16}},da2)
        else
            return da2
        end
    elseif eltype_old <: AbstractFloat
        # first test if the floats are integer numbers
        if sum(isinteger.(collect(skipmissing(da2))) .== false) == 0
            return acompress(convert(Vector{Union{Missing,Int64}},da2))
        end

        # get minimum and maximum values
        if eltype_old == Float64 && minimum(collect(skipmissing(da2))) >= floatmin(Float32) && maximum(collect(skipmissing(da2))) <= floatmax(Float32)
            return convert(Vector{Union{Missing,Float32}},da2)
        else
            return da2
        end
    end
end


function atype(df::DataFrame,v::Symbol)
    # Array type = DA for DataArray, CA for Categorical Array, and UV for Union Vector
    if isdefined(Main,:CategoricalArrays) && typeof(df[!,v]) <: CategoricalArray
        return string("Categorical (", replace(string(eltype(df[!,v].refs)),"UInt" => ""), ")")
    elseif isdefined(Main,:DataArrays) && typeof(df[!,v]) <: DataArray
         return "DataArray"
    elseif isdefined(Main,:PooledArrays) && typeof(df[!,v]) <: PooledArray
         return "PooledArray"
    elseif isa(eltype(df[!,v]),Union)
        return "Union Vector" # Union Vector
    else
        return "Vector"
    end
end

function etype(df::DataFrame,v::Symbol)
    # Eltype
    if typeof(df[!,v]) <: CategoricalArray
        eltyp = string(eltype(df[!,v].pool.index))
        if in(eltyp,["String","AbstractString"])
            eltyp = string("Str",getmaxwidth(df[!,v].pool.index))
        end
    else
        eltyp = string(nonmissingtype(eltype(df[!,v])))
        if in(eltyp,["String","AbstractString"])
            eltyp = string("Str",getmaxwidth(df[!,v]))
        elseif eltyp == "Dates.Date"
            eltyp = "Date"
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
    desc(df::DataFrame,varnames::Symbol...; labels::Union{Nothing,Label}=nothing)

Displays variables in a dataframe much like `showcols`. It can display additional
attributes such as variable labels, value labels and display formats (not used in Julia)
if an optional `labels` is specified. It mimics Stata's `describe` command.
`labels` is automatically converted from a stata file by `read_stata` function. Or one can
be easily created as described in [Labels](https://github.com/mwsohn/Labels.jl).
"""
function desc(df::DataFrame,varnames::Symbol...; labels::Union{Nothing,Label}=nothing)

    if length(varnames) == 0
        varnames = names(df)
    end

    varlen = zeros(Int,size(df,2)) # length of variable names
    lablen = zeros(Int,size(df,2)) # length of value labels
    for (i,v) in enumerate(varnames)
        varlen[i] = length(string(v))
        lablen[i] = labels != nothing && lblname(labels,v) != nothing ? length(string(lblname(labels,v))) : 0
    end

    # width for variable names
    maxval = maximum(varlen)
    maxval = maxval < 8 ? 8 : maxval

    # width for variable types
    maxatype = 16
    maxeltype = 7
    maxmiss = 7

    # width for value label names
    maxlab = maximum(lablen)
    maxlab = maxlab < 11 ? 11 : maxlab

    # width for display formats
    # maxformat = maximum(forlen)
    # maxformat = maxformat < 6 ? 6 : maxformat

    # number of variables
    numvar = length(varnames)

    # width for the variable index - minimum 3 spaces
    maxobs = length(string(numvar))
    maxobs = maxobs < 3 ? 3 : maxobs
    print(lpad("Num",maxobs),"  ")

    # variable name
    print(rpad("Variable",maxval),"  ")

    # type of array
    print(rpad("Atype",maxatype),"  ")

    # type (eltype)
    print(rpad("Eltype",maxeltype),"  ")

    # percent missing
    print(rpad("Missing",maxmiss),"  ")

    if labels != nothing
        # label
        print("Description")
    end

    print("\n")

    # dashes for a line by itself -- assume 30 characters for "label"
    numdashes = maxobs+maxval+maxatype+maxmiss+maxeltype+4
    if labels == nothing
        println(repeat("-",numdashes+4))
    else
        println(repeat("-",numdashes+30))
    end

    nrows,ncols = size(df)

    # output dataframe
    dfv = DataFrame(Variable = varnames)
    dfv[!,:ArrayType] = Vector{String}(undef,ncols)
    dfv[!,:Eltype] = Vector{String}(undef,ncols)
    dfv[!,:Missing] = Vector{String}(undef,ncols)
    dfv[!,:Lblname] = Vector{String}(undef,ncols)
    dfv[!,:Description] = Vector{String}(undef,ncols)


    for (i,v) in enumerate(varnames)

        # variable name
        varstr = string(v)

        # Array Type
        dfv[i,:ArrayType] = atype(df,v)

        # Eltype
        dfv[i,:Eltype] = etype(df,v)

        # percent missing
        nmiss = sum(Missings.ismissing.(df[!,v]))
        dfv[i,:Missing] = string(round(100 * nmiss/nrows,digits=1),"%")

        print(lpad(string(i),maxobs),"  ",rpad(varstr,maxval),"  ",
            rpad(dfv[i,:ArrayType],maxatype),"  ",
            rpad(dfv[i,:Eltype],maxeltype),"  ",
            lpad(dfv[i,:Missing],maxmiss),"  ")

        if labels != nothing
            dfv[i,:Lblname] = lblname(labels,v) == nothing ? "" : string(lblname(labels,v))
            dfv[i,:Description] = varlab(labels,v)
            print(dfv[i,:Description])
        end

        print("\n")
    end

    return dfv
end

function getmaxwidth(s::AbstractArray)
    return  maximum(length.(collect(skipmissing(s))))
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

    dslist = Array{Symbol,1}()

    for v in names(df)
        if typ == String && eltype(df[v]) == String
            if length(args) == 0
                push!(dslist,v)
                continue
            end

            maxlen = getmaxwidth(df[v])
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

    for v in names(df)
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
    df[:_____obs_____] = collect(1:size(df,1))
    done = zeros(Int8,size(df,1))
    for subdf in groupby(df, groupvars)
        done[subdf[1,:_____obs_____]]=1
    end
    deletecols!(df,:_____obs_____)
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
function duplicates(df::DataFrame, args::Symbol... ; cmd::Symbol = :report)

    if in(cmd,[:report,:drop,:tag]) == false
        error("`cmd = `", cmd, "` is not supported.")
    end

    # if args are not specified, use all variables
    nargs = length(args)
    if nargs == 0
        args = tuple(names(df)...)
    end

    dfx = by(df,collect(args),x->size(x,1))

    rename!(dfx,:x1 => :__dups)
    if cmd == :report
        na = freqtable(dfx, :__dups)
        setdimnames!(na,"copies",1)
        return na
    end

    # substract 1 from :x in dfx (we are reporting 0 for unique observations)
    dfx[:__dups] .-= 1
    df = join(df, dfx, on = collect(args), kind = :left)

    if cmd == :tag
        return df[:__dups]
    end

    # if cmd == :drop
    ba = [x == 1 ? true : false for x in pickone(df,collect(args))]
    return df[ba,collect(names(df[1:end-1]))]
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
function dfmerge(df1::DataFrame,df2::DataFrame,linkers::Union{Symbol,Vector};kind::Symbol = :outer)
    if in(:_merge,names(df1))
        error("`:_merge' exists in the first dataframe")
    end
    if in(:_merge,names(df2))
        error("`:_merge' exists in the second dataframe")
    end

    df1[:___mergeleft___] = ones(Int8,size(df1,1))
    df2[:___mergeright___] = ones(Int8,size(df2,1))

    df_merged = join(df1,df2,on = linkers,kind=kind)
    df_merged[:_merge] = zeros(Int8,size(df_merged,1))
    for i = 1:size(df_merged,1)
        if ismissing(df_merged[i,:___mergeright___])
            df_merged[i,:_merge] = 1
        elseif ismissing(df_merged[i,:___mergeleft___])
            df_merged[i,:_merge] = 2
        elseif df_merged[i,:___mergeleft___] == 1 && df_merged[i,:___mergeright___] == 1
            df_merged[i,:_merge] = 3
        end
    end
    delete!(df_merged,[:___mergeleft___,:___mergeright___])
    return df_merged
end

"""
    dfappend(df1::DataFrame,df2::DataFrame)

Produces a dataframe with two sources stacked together. It is simply a wrapper
for `[df1;df2]` operation to micmic `append using ...` command in Stata.
"""
function dfappend(df1::DataFrame,df2::DataFrame)
    return [df1;df2]
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
            rename!(df,nm => Symbol(newname))
        end
    end
end

vidx(df::DataFrame,varname::Symbol) = findfirst(x->x == varname, names(df))

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
    da2 = Vector{Union{Missing,T}}(length(da))

    for i in 1:length(da)
        da2[i] = ismissing(da[i]) ? missing : parse(T,da[i])
    end

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

# """
#     rowsum(df::DataFrame)
#
# Creates a DataArray that contains the row total of all values on the same row of `df`.
# If one of the columns contain an NA value on a row, an NA value will be returned for that
# row. This function emulates Stata's `egen rowsum = rowtotal(var1 - var3)`.
#
# ```
# julia>df[:rowsum] = df[[:var1,:var2,:var3]]
# ```
#
# If the position numbers for `:var1` (e.g., 4), `:var2` (5), `:var3` (6) are known and consecutive,
# you can specify them as follows:
#
# ```
# julia>df[:rowsum] = df[collect(4:6)]
# ```
# """
# function rowsum(df::DataFrame)
#
#     isfloat = false
#     for i in 1:size(df,2)
#         if eltype(df[i]) <: AbstractFloat
#             isfloat = true
#         end
#     end
#
#     if isfloat
#         da = zeros(Union{Missing,Float64},size(df,1))
#     else
#         da = zeros(Union{Missing,Int64},size(df,1))
#     end
#
#     ba = completecases(df)
#     for i = 1:size(df,1)
#         if ba[i] == false
#             da[i] = missing
#         else
#             for j = 1:size(df,2)
#                 da[i] += df[i,j]
#             end
#         end
#     end
#     return da
# end

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
    rapply(df::DataFrame,r::OrderedDict)

Creates a recoded vector of values according to the rule set specified as an ordered dictionary.
Keys in the rule set must be the recoded value and the values must be the rules.
The rule must be an expression or a string that return a boolean value (`true` or `false`).

## Example

julia> ruleset = OrderedDict(
            1 => :( df[:race] .== 1 && df[:hispanic] .== 0),
            2 => :( df[:race] .== 2 && df[:hispanic] .== 0),
            3 => :( df[:hispanic] .== 1),
            4 => :( df[:hispanic] .== 0 && in.(df[:race],[1,2]) == false)
       )
OrderedDict{Int64,Expr} with 4 entries:
  1 => :(df[:race] .== 1 && df[:hispanic] .== 0)
  2 => :(df[:race] .== 2 && df[:hispanic] .== 0)
  3 => :(df[:hispanic] .== 1)
  4 => :(df[:hispanic] .== 0 && in.(df[:race], [1, 2]) == false)

julia> recode(adf,ruleset)


"""
function rapply(df::DataFrame,r::OrderedDict)

    # values
    vals = sort(collect(keys(r)))

    # types
    vtyp = eltype(vals)

    # empty Vector
    vec = Vector{Union{vtyp,Missing}}(undef,size(df,1))

    # go through the rules and assign values
    for v in vals

        if typeof(r[v]) == String
            ba = eval(parse(r[v]))
        elseif typeof(r[v]) == Expr
            ba = eval(r[v])
        else
            error(typeof(r[v]), " is not an allowed type for rule ", string(r[v]))
        end

        for i in findall(x->x==true,ba)
            vec[i] = v
        end
    end

    return vec
end
