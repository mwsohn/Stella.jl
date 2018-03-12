function saved_memory(dt,n)

    # total size and levels
    tsize = 0
    s = Set()
    for i=1:dt.rows
        tsize += length(dt.data[n][i]).value
        if dt.data[n][i].hasvalue
            push!(s, dt.data[n][i].value)
        end
    end

    # for refs size, assume UInt32
    refsize = 4

    rsize = sum([length(x) for x in s]) + refsize*dt.rows

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
            df[dt.headers[i]] = Union{Missing,Int32}[x.hasvalue ? dateoffset + Dates.Day(x.value) : missing for x in dt.data[i]]
        elseif dt.types[i] <: Integer && dt.formats[i] in ("%tc","%tC")
            df[dt.headers[i]] = Union{Missing,Int64}[x.hasvalue ? datetimeoffset + Dates.Millisecond(x.value) : missing for x in dt.data[i]]
        elseif dt.types[i] <: AbstractString
            if saved_memory(dt,i)
                df[dt.headers[i]] = categorical(String[x.hasvalue ? x.value : "" for x in dt.data[i]],true)
            else
                df[dt.headers[i]] = String[x.hasvalue ? x.value : "" for x in dt.data[i]]
            end
        elseif dt.types[i] == Int8 && dt.val_label_keys[i] != ""
            # variable is Int8 type and has label name, convert it to categorical array
            df[dt.headers[i]] = categorical(Union{Missing,dt.types[i]}[x.hasvalue ? x.value : missing for x in dt.data[i]],true)
        else
            df[dt.headers[i]] = Union{Missing,dt.types[i]}[x.hasvalue ? x.value : missing for x in dt.data[i]]
        end
    end
    return df
end
import DataFrames.DataFrame
DataFrame(x::ReadStat.ReadStatDataFrame) = readstat2dataframe(x)

function read_stata(fn::String,verbose=false)

    dt = ReadStat.read_dta(fn)
    df = readstat2dataframe(dt,verbose)

    # labels
    varlab = Dict()
    lblname = Dict()
    formatlab = Dict()

    for i=1:dt.columns
        varlab[dt.headers[i]] = dt.labels[i]
        lblname[dt.headers[i]] = dt.val_label_keys[i]
        formatlab[dt.headers[i]] = dt.formats[i]
    end

    vallab = Dict()
    for k in collect(keys(dt.val_label_dict))
        vallab[k] = dt.val_label_dict[k]
    end

    label = Dict()
    label["variable"] = varlab
    label["label"] = lblname
    label["formats"] = formatlab
    label["values"] = vallab

    return df,label
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
        if length(df[v]) == sum(df[v].na)
            delete!(df,v)
            println(v, " was empty, now deleted")
            continue
        end

        # get the original eltype
        eltype_old = eltype(df[v])

        df[v] = acompress(df[v])

        if eltype_old != eltype(df[v])
            println(v, " was ", eltype_old, ", now ", eltype(df[v]))
        end
    end
end

function acompress(da::AbstractVector)
    # get the original eltype
    eltype_old = eltype(da)

    # if the original eltype is a Boolean or a String
    # then do not perform any compression
    # if eltype_old == Bool
    #     warn(da," is a Boolean. No need to compress.")
    # elseif eltype_old <: AbstractString
    #     warn(da," is a String. No need to compress.")
    # elseif length(dropna(da)) == 0
    #     warn(da," is an empty column. No need to compress")
    # end

    if  eltype_old <: Integer
        # get minimum and maximum values
        varmin = minimum(dropna(da))
        varmax = maximum(dropna(da))
        # test if Int8
        if eltype_old != Int8 && varmin >= typemin(Int8) && varmax <= typemax(Int8)
            return convert(Vector{Union{Missing,Int8}},da)
        elseif eltype_old != Int16 && varmin >= typemin(Int16) && varmax <= typemax(Int16)
            return convert(Vector{Union{Missing,Int16}},da)
        else
            return da
        end
    elseif eltype_old <: AbstractFloat
        # first test if the floats are integer numbers
        if sum(isinteger.(dropna(da)).==false) == 0
            return compress(convert(Vector{Union{Missing,Int64}},da))
        end

        # get minimum and maximum values
        if eltype_old != Float32 && minimum(dropna(da)) >= typemin(Float32) && maximum(dropna(da)) <= typemin(Float32)
            return convert(Vector{Union{Missing,Float32}},da)
        else
            return da
        end
    end
end

"""
    desc(df::DataFrame,varnames::Symbol...; label_dict::Union{Void,Dict}=nothing)

Displays variables in a dataframe much like `showcols`. It can display additional
attributes such as variable labels, value labels and display formats (not used in Julia)
if an optional `label_dict` is specified. It mimics Stata's `describe` command.
`label_dict` is automatically converted from a stata file by `read_stata` function. Or one can
be easily created as follows:

## Label Dictionary Format

```jldoctest
label = Dict()
label["variable"] = Dict(
    "id" => "Identification Number",
    "age" => "Age in year at time of interview"),
    "sex" => "Respondent sex"
    )

label["value"] = Dict()
label["value"]["sexlabel"] = Dict(
    1 => "Female",
    2 => "Male"
    )

label["label"] = Dict(
    "sex" => "sexlabel"
    )
```
"""
function desc(df::DataFrame,varnames::Symbol...;label_dict::Union{Void,Dict}=nothing)

    if length(varnames) == 0
        varnames = names(df)
    end

    varlab = label_dict != nothing && haskey(label_dict,"variable") ? label_dict["variable"] : Dict()
    lablab = label_dict != nothing && haskey(label_dict,"label") ? label_dict["label"] : Dict()
    forlab = label_dict != nothing && haskey(label_dict,"format") ? label_dict["format"] : Dict()

    varlen = zeros(Int,size(df,2)) # length of variable names
    lablen = zeros(Int,size(df,2)) # length of value labels
    forlen = zeros(Int,size(df,2)) # length of display formats
    for (i,v) in enumerate(varnames)
        varstr = string(v)
        varlen[i] = length(varstr)
        if label_dict != nothing
            if haskey(lablab,varstr)
                tmplen = length(lablab[varstr])
                lablen[i] = length(lablab[varstr])
            end
            if haskey(forlab,v)
                tmplen = length(forlab[v])
                forlen[i] = length(forlab[v])
            end
        end
    end

    # width for variable names
    maxval = maximum(varlen)
    maxval = maxval < 8 ? 8 : maxval

    # width for variable types
    maxatype = 5
    maxeltype = 7
    maxmiss = 7

    # width for value label names
    maxlab = maximum(lablen)
    maxlab = maxlab < 11 ? 11 : maxlab

    # width for display formats
    maxformat = maximum(forlen)
    maxformat = maxformat < 6 ? 6 : maxformat

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

    if label_dict != nothing
        # # value label
        # print(rpad("Value Label",maxlab),"  ")
        #
        # # format
        # print(rpad("Format",maxformat),"  ")

        # label
        print("Label")
    end

    print("\n")

    # dashes for a line by itself -- assume 30 characters for "label"
    numdashes = maxobs+maxval+maxatype+maxmiss+maxeltype+4
    if label_dict == nothing
        println(repeat("-",numdashes+4))
    else
        println(repeat("-",numdashes+30))
    end

    nrows = size(df,1)
    for (i,v) in enumerate(varnames)

        # variable name
        varstr = string(v)

        # Array type = DA for DataArray, CA for Categorical Array, and UV for Union Vector
        if isa(eltype(df[v]),Union)
            atyp = "UV" # Union Vector
        elseif typeof(df[v]) <: CategoricalArray
            atyp = "CA"
        else
            atyp = "VC"
        end

        # Eltype
        if typ <: CategoricalArray
            eltyp = UInt32
        else
            eltyp = string(Missings.T(eltype(df[v])))
        end
        # eltyp = string(Missings.T(eltype(df[v])))

        if in(eltyp,["String","AbstractString"])
            eltyp = string("Str",getmaxlength(df[v]))
        end

        # percent missing
        nmiss = typeof(df[v]) <: DataArray ? sum(isna.(df[v])) : sum(Missings.ismissing.(df[v]))
        pmiss = string(round(100 * nmiss/nrows,1),"%")

        print(lpad(string(i),maxobs),"  ",rpad(varstr,maxval),"  ",lpad(atyp,maxatype),"  ",rpad(eltyp,maxeltype),"  ",lpad(pmiss,maxmiss),"  ")

        if label_dict != nothing
            # if haskey(lablab,v)
            #     print(rpad(lablab[v],maxlab),"  ")
            # else
            #     print(repeat(" ",maxlab+2))
            # end
            #
            # if haskey(forlab,v)
            #     print(rpad(forlab[v],maxformat),"  ")
            # else
            #     print(repeat(" ",maxformat+2))
            # end
            #
            if haskey(varlab,v)
                print(varlab[v])
            end
        end

        print("\n")
    end
end

function getmaxlength(s::AbstractArray)

    if length(collect(skipmissing(s))) == 0
        return 0
    end

    # function does not work on empty arrays
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
```jldoctest
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

            maxlen = getmaxlength(df[v])
            if length(args) == 1
                if args[1] <= maxlen
                    push!(dslist,v)
                end
            elseif length(args) == 2
                if args[1] <= maxlen <= args[2]
                    push!(dslist,v)
                end
            end
        elseif eltype(df[v]) == typ
            push!(dslist,v)
        elseif (typ == Integer && eltype(df[v]) <: Integer) || (typ == AbstractFloat && eltype(df[v]) <: AbstractFloat) || (typ == Number && eltype(df[v]) <: Number)
            push!(dslist,v)
        end
    end

    return dslist

end
function ds(df::DataFrame,re::Regex)
    dslist = Array{Symbol,1}()

    for v in names(df)
        if ismatch(re,string(v))
            push!(dslist,v)
        end
    end

    return dslist
end

"""
    duplicates(df::DataFrame, args::Symbol...; cmd::Symbol = :report)

Reports, tags, or deletes duplicates in a DataFrame for one or more columns. Use `cmd` to
request one of the actions: `:report` for getting a frequency table for the number of duplicates;
`:tag` for identifying rows with duplicate values; and `:drop` for dropping all rows with duplicate values
except for the first row.
"""
function duplicates(df::DataFrame, args::Symbol... ; cmd::Symbol = :report) #; exclude::Array{Union{Symbol,Int},1} = [], include::Array{Union{Symbol,Int},1} = [] )

    if in(cmd,[:report,:drop,:tag]) == false
        error("`cmd = `", cmd, "` is not supported.")
    end

    nargs = length(args)
    if nargs == 0
        args = tuple(names(df)...)
    end

    dfx = by(df,collect(args),x->size(x,1))

    rename!(dfx,:x1,:__dups)
    if cmd == :report
        na = freqtable(dfx, :__dups)
        setdimnames!(na,"copies",1)
        return na
    end

    # substract 1 from :x in dfx (we are reporting 0 for unique observations)
    dfx[:__dups] -= 1
    df = join(df, dfx, on = collect(args), kind = :left)

    if cmd == :tag
        return df[:__dups]
    end

    # cmd == :drop
    ba = [x == 1 ? true : false for x in pickone(df,collect(args))]
    return df[ba,collect(names(df[1:end-1]))]
end

"""
    dfsample(df::AbstractDataFrame,num::Union{Int64,Float64})

Creates a dataframe with a randomly selected sample from the input dataframe. `num`
specifies the sample size. If `num` is an integer, it indicates the number of rows to be selected.
If it is a float, it indicates the percentage of the input dataframe to be randomly selected.
Use `srand()` to set a seed before running `dfsample()`.

##Example
To select 100 rows, use

```jldoctest
julia> df2 = dfsample(df,100)
```

To select 20% of the original dataframe, use

```jldoctest
julia> df2 = dfsample(df,20.)
```

or

```jldoctest
julia> df2 = dfsample(df,.2)
```

"""
function dfsample(df::AbstractDataFrame,num::Union{Int64,Float64})

    df2 = DataFrame()
    df2[:___ran_num___] = rand(Uniform(),size(df,1))
    df2[:___obs___] = collect(1:size(df,1))

    if typeof(num) <: AbstractFloat

        if 1. <= num < 100.
            num /= 100.
        elseif num >= 100
            error("A percentage value in floating point number (e.g., 10.0 or .1) is required.")
        end

        num2 = round(Int,num*size(df,1))
    else
        # number of obserations
        if num > size(df,1)
            error(num, " is greater than the total number of rows in the input dataframe.")
        end
        num2 = num
    end

    sort!(df2,cols=[:___ran_num___])
    a = convert(Vector,df2[1:num2,:___obs___])
    return df[sort(a),:]
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
        if isna(df_merged[i,:___mergeright___])
            df_merged[i,:_merge] = 1
        elseif isna(df_merged[i,:___mergeleft___])
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
    classify(da::DataArray,thresholds::Vector; lower::Bool = false)

Produces a categorical variable based on a data array and a vector of thresholds.
Use `lower = true` to classify the threshold value to the lower class.
"""
function classify(da::AbstractVector,thresholds::Vector; lower::Bool = false)
    da2 = Vector{Union{Missing,Int8}}(length(da))
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
    addvars(fm::Formula,v::Vector{Symbol})

Adds covariates to an existing formula. `fm` is an object of `Formula` type
created by `@formula` macro. `v` is a vector of Symbols. This function is useful to
build hierarchically nested models.

### Example

```jldoctest
julia> using DataFrames

julia> fm = @formula(income ~ age + male)
Formula: income ~ age + male

julia> fm2 = addvars(fm,[:race, :educ])
Formula: income ~ age + male + race + educ
```

"""
function addvars(fmm::Formula,v::Vector{Symbol})

    # create a new Formula object
    fm = deepcopy(fmm)

    # if fm.rhs is a Symbol, convert it to an Expr first
    if typeof(fm.rhs) == Symbol
        tmpvar = fm.rhs
        fm.rhs = :()
        push!(fm.rhs.args,:+)
        push!(fm.rhs.args,tmpvar)
        fm.rhs.head = :call
    end

    for i=1:length(v)
        push!(fm.rhs.args,v[i])
    end

    return fm
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
            rename!(df,nm,Symbol(newname))
        end
    end
end

vidx(df::DataFrame,varname::Symbol) = findfirst(varname, names(df))

"""
    destring(da::DataArray;force=true)
    destring(df::DataFrmae,strvar::Symbol;force=true)
    destring!(df::DataFrame,strvars; newvars::Vector{Symbol} = [], force=true, replace=false)

Convert a string DataArray to a numeric DataArray. Use `force = true` to coerce conversion of alphanumeric strings to
`NA` values. If `force = false`, any DataArray containing non-numeric values are not converted.
Use `replace = true` in `destring!` to replace the original string DataArray with a new converted numeric DataArray.
If `replace` option is specified, `newvars` array is ignored.
"""
function destring(da::AbstractArray; force=true)
    if length(da) == 0
        error(da,"Input data array is empty!")
    end
    if eltype(da) <: Number
        error(da," is a numeric DataArray.")
    end

    # check if the values include any alphabetic characters or decimals
    isfloat = false
    alpha = false
    da_safe = dropna(da)
    for i in length(da_safe)
        if sum([isalpha(x) for x in da_safe[i]]) > 0
            alpha = true
        end
        if ismatch(r"[,0-9]*\.[0-9]+",da_safe[i])
            isfloat = true
        end
    end

    if alpha && force == false
        error(arg," contains alphabetic letters. Use 'force=true' option to coerce conversion.")
    end

    T = isfloat ? Float64 : Int64
    da2 = Vector{Union{Missing,T}}(length(da))

    for i in 1:length(da)
        da2[i] = ismissing(da[i]) ? missing : parse(T,da[i])
    end

    return compress(da2)
end
destring(df::DataFrame,strvar::Symbol; force=true) = destring(df[strvar],force=force)
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
    rowsum(df::DataFrame)

Creates a DataArray that contains the row total of all values on the same row of `df`.
If one of the columns contain an NA value on a row, an NA value will be returned for that
row. This function emulates Stata's `egen rowsum = rowtotal(var1 - var3)`.

```
julia>df[:rowsum] = df[[:var1,:var2,:var3]]
```

If the position numbers for `:var1` (e.g., 4), `:var2` (5), `:var3` (6) are known and consecutive,
you can specify them as follows:

```
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
        da = zeros(Union{Missing,Float64},size(df,1))
    else
        da = zeros(Union{Missing,Int64},size(df,1))
    end

    ba = completecases(df)
    for i = 1:size(df,1)
        if ba[i] == false
            da[i] = missing
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

```
julia>df[:rowmean] = rowstat(df[[:var1,:var2,:var3]],mean)
```

If the position numbers for `:var1` (e.g., 4), `:var2` (5), `:var3` (6) are known and consecutive,
you can specify them as follows:

```
julia>df[:rowstd] = rowstat(df[collect(4:6)],std)
```
"""
function rowstat(df::DataFrame,func::Function)

    da = zeros(Float64,size(df,1))
    ta = Array{Float64,1}(size(df,2))

    for i = 1:size(df,1)
        k=0
        for j = 1:size(df,2)
            if isna(df[i,j]) == false
                k += 1
                ta[k] = df[i,j]
            end
        end
        if k == 0
            da[i] = missing
        else
            tmpfloat = func(ta[1:k])
            da[i] = isnan(tmpfloat) ? missing : tmpfloat
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

```
julia> df[:agecat] = xtile(df[:age], nq = 3)
```
"""
function xtile(da::AbstractArray ; nq::Int = 4, cutoffs::Union{Void,AbstractVector} = nothing)

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
	    cutoffs = nquantile(dropmissing(da),nq)
    elseif length(cutoffs) == nq - 1
        cutoffs = vcat(minimum(dropmissing(da)), cutoffs, maximum(dropmissing(da)))
    else
        error("`cutoffs` vector length is not consistent with `nq`. It must be 1 greater or 1 less than `nq`.")
    end

	return [ismissing(x) ? missing : qval(x,cutoffs) for x in da]
end
xtile(df::DataFrame,arg::Symbol; nq::Int = 4, cutoffs::Union{Void,AbstractVector} = nothing) = xtile(df[arg], nq = nq, cutoffs = cutoffs)
