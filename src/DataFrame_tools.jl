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
        if length(skipmissing(df[v])) == 0
            delete!(df,v)
            println(v, " was empty, now deleted")
            continue
        end

        # get the original eltype
        eltype_old = eltype(df[v])

        df[v] = compress(df[v])

        if eltype_old != eltype(df[v])
            println(v, " was ", eltype_old, ", now ", eltype(df[v]))
        end
    end
end

function compress(da::AbstractVector)
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
        varmin = minimum(skipmissing(da))
        varmax = maximum(skipmissing(da))
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
        if sum(isinteger.(skipmissing(da)).==false) == 0
            return compress(convert(Vector{Union{Missing,Int64}},da))
        end

        # get minimum and maximum values
        if eltype_old != Float32 && minimum(skipmissing(da)) >= typemin(Float32) && maximum(skipmissing(da)) <= typemin(Float32)
            return convert(Vector{Union{Missing,Float32}},da)
        else
            return da
        end
    end
end

function get_disp_length(dat::AbstractVector; precision = 2)

    vtype = eltype(dat)
    if vtype == Date
        return 10 # 10/10/2017
    end

    if vtype == DateTime # 10/10/2017 14:20
        return 16
    end

    if vtype <: Number
        da = skipmissing(dat)
        _max = size(da,1) == 0 ? 0 : maximum(da)
        for k = 1:30
            if _max < 10^k
                if vtype <: Integer || precision == 0
                    return k
                else
                    return k + precision + 1
                end
            end
        end
        error("Cannot display numbers that are 10³¹ long.")
    end

    if vtype <: AbstractString
        return Stella.getmaxlength(dat)
    end

    error(eltype(dat)," is not recognized.")
end


function dflist(df::DataFrame; precision = 2, n = 10)

    nvars = size(df,2)
    vars = names(df)

    lenvec = [ get_disp_length(df[y], precision = precision) for y in vars ]
    # println(lenvec)

    # change lenvec to account for column names up to 8 characters
    vtypes =eltypes(df)

    #headers
    for i = 1:nvars
        if vtypes[i] <: Number || vtypes[i] == Date || vtypes[i] == DateTime
            print(prepend_spaces(string(vars[i]),lenvec[i]))
        else
            print(append_spaces(string(vars[i]),lenvec[i]))
        end

        if i < nvars
            print("  ")
        end
    end
    print("\n")

    if n==0
        n = size(df,1)
    end

    for ln = 1:n
        for i = 1:nvars
            if isna(df[ln,i])
                if vtypes[i] <: AbstractString
                    print(append_spaces("NA",lenvec[i]))
                else
                    print(prepend_spaces("NA",lenvec[i]))
                end
                if ln < nvars
                    print("  ")
                end
                continue
            end
            if vtypes[i] == DateTime
                print(Dates.format(df[ln,i],"mm/dd/yyyy HH:MM"))
            elseif vtypes[i] == Date
                print(Dates.format(df[ln,i],"mm/dd/yyyy"))
            elseif vtypes[i] <: AbstractString
                print(append_spaces(df[ln,i],lenvec[i]))
            elseif vtypes[i] <: Integer
                print(prepend_spaces(string(df[ln,i]),lenvec[i]))
            elseif vtypes[i] <: AbstractFloat
                if precision == 0
                    nstr = @sprintf("%.0f",df[ln,i])
                elseif precision == 1
                    nstr = @sprintf("%.1f",df[ln,i])
                elseif precision == 2
                    nstr = @sprintf("%.2f",df[ln,i])
                elseif precision == 3
                    nstr = @sprintf("%.3f",df[ln,i])
                elseif precision == 4
                    nstr = @sprintf("%.4f",df[ln,i])
                elseif precision == 5
                    nstr = @sprintf("%.5f",df[ln,i])
                elseif precision == 6
                    nstr = @sprintf("%.6f",df[ln,i])
                elseif precision == 7
                    nstr = @sprintf("%.7f",df[ln,i])
                elseif precision == 8
                    nstr = @sprintf("%.8f",df[ln,i])
                else
                    error("Can't display that much precision.")
                end
                print(prepend_spaces(nstr,lenvec[i]))
            end
            if i < nvars
                print("  ")
            end
        end
        print("\n")
    end
end


"""
    desc(df::DataFrame; label_dict::Union{Void,Dict}=nothing)

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
function desc(df::DataFrame,varnames::Vector=[];label_dict::Union{Void,Dict}=nothing)

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
    maxtype = 7

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
    print(prepend_spaces("Num",maxobs),"  ")

    # variable name
    print(append_spaces("Variable",maxval),"  ")

    # type (eltype)
    print(append_spaces("Type",maxtype),"  ")

    if label_dict != nothing
        # value label
        print(append_spaces("Value Label",maxlab),"  ")

        # format
        print(append_spaces("Format",maxformat),"  ")

        # label
        print("Label")
    end

    print("\n")

    # dashes for a line by itself -- assume 20 characters for "label"
    if label_dict == nothing
        println(repeat("-",maxobs+maxval+maxtype+4))
    else
        println(repeat("-",maxobs+maxval+maxtype+maxlab+maxformat+30))
    end

    for (i,v) in enumerate(varnames)

        eltyp = string(eltype(df[v]))

        if in(eltyp,["String","AbstractString"])
            eltyp = string("Str",getmaxlength(df[v]))
        end

        varstr = string(v)
        print(prepend_spaces(string(i),maxobs),"  ",append_spaces(varstr,maxval),"  ",append_spaces(eltyp,maxtype),"  ")

        if label_dict != nothing
            if haskey(lablab,v)
                print(append_spaces(lablab[v],maxlab),"  ")
            else
                print(repeat(" ",maxlab+2))
            end

            if haskey(forlab,v)
                print(append_spaces(forlab[v],maxformat),"  ")
            else
                print(repeat(" ",maxformat+2))
            end

            if haskey(varlab,v)
                print(varlab[v])
            end
        end

        print("\n")
    end
end

function getmaxlength(s::DataArray)
    if sum(s.na) < size(s,1)
        return maximum(length.(dropna(s)))
    end
    return 0
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
    dfsample(df::DataFrame,num::Union{Int64,Float64})

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
function dfsample(df::DataFrame,num::Union{Int64,Float64})

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

Produces a categorical variable based on a data array and a vector of throsholds.
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
