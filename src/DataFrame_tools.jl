#############################################################################
#
# DataFrame tools
#
#############################################################################
"""
    dfcompress!(df::DataFrame)

Reduce `df`'s memory use by changing the eltype of each DataArray in the `df` to
the type that can accommodate the larget and the smallest values within the same
integer or float class.
"""
function dfcompress!(df::DataFrame)
    for v in names(df)

        # if DataArray is empty after all NAs are dropped
        # drop it from the df
        if length(dropna(df[v])) == 0
            delete!(df,v)
            println(v, " was empty, now deleted")
            continue
        end

        # get the original eltype
        eltype_old = eltype(df[v])

        df[v] = dacompress(df[v])

        if eltype_old != eltype(df[v])
            println(v, " was ", eltype_old, ", now ", eltype(df[v]))
        end
    end
end

function dacompress(da::DataArray)
    # get the original eltype
    eltype_old = eltype(da)

    # if the original eltype is a Boolean or a String
    # then do not perform any compression
    if eltype_old == Bool
        error(da," is a Boolean. No need to compress.")
    elseif eltype_old <: AbstractString
        error(da," is a String. No need to compress.")
    elseif length(dropna(da)) == 0
        error(da," is an empty column. No need to compress")
    end

    if  eltype_old <: Integer
        # get minimum and maximum values
        varmin = minimum(dropna(da))
        varmax = maximum(dropna(da))
        # test if Int8
        if eltype_old != Int8 && varmin >= typemin(Int8) && varmax <= typemax(Int8)
            return convert(DataArray{Int8,1},da)
        elseif eltype_old != Int16 && varmin >= typemin(Int16) && varmax <= typemax(Int16)
            return convert(DataArray{Int16,1},da)
        else
            return da
        end
    elseif eltype_old <: AbstractFloat
        # get minimum and maximum values
        varfmin = minimum(dropna(da))
        varfmax = maximum(dropna(da))
        if eltype_old != Float32 && varfmin >= typemin(Float32) && varfmax <= typemin(Float32)
            return convert(DataArray{Float32,1},da)
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
        _max = maximum(dropna(dat))
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
        return getmaxlength(dat)
    end

    error(eltype(dat)," is not recognized.")
end


function dflist(df::DataFrame,args::Symbol...; precision = 2, n = 10)

    nargs = length(args)

    lenvec = [ get_disp_length(df[y], precision = precision) for y in 1:nargs ]
    println(lenvec)

    # change lenvec to account for column names up to 8 characters

    vtypes =eltypes(df[collect(args)])
    println(vtypes)

    #headers
    for i = 1:nargs
        if vtypes[i] <: Number || vtypes[i] == Date || vtypes[i] == DateTime
            print(prepend_spaces(string(args[i]),lenvec[i]))
        else
            print(append_spaces(string(args[i]),lenvec[i]))
        end

        if i < nargs
            print("  ")
        end
    end
    print("\n")

    if n==0
        n = size(df,1)
    end

    for ln = 1:n
        for i = 1:nargs
            if isna(df[ln,args[i]])
                if vtypes[i] <: AbstractString
                    print(append_spaces("NA",lenvec[i]))
                else
                    print(prepend_spaces("NA",lenvec[i]))
                end
                if ln < nargs
                    print("  ")
                end
                continue
            end
            if vtypes[i] == DateTime
                print(Dates.format(df[ln,args[i]],"mm/dd/yyyy HH:MM"))
            elseif vtypes[i] == Date
                print(Dates.format(df[ln,args[i]],"mm/dd/yyyy"))
            elseif vtypes[i] <: AbstractString
                print(append_spaces(df[ln,args[i]],lenvec[i]))
            elseif vtypes[i] <: Integer
                print(prepend_spaces(string(df[ln,args[i]]),lenvec[i]))
            elseif vtypes[i] <: AbstractFloat
                if precision == 0
                    nstr = @sprintf("%.0f",df[ln,args[i]])
                elseif precision == 1
                    nstr = @sprintf("%.1f",df[ln,args[i]])
                elseif precision == 2
                    nstr = @sprintf("%.2f",df[ln,args[i]])
                elseif precision == 3
                    nstr = @sprintf("%.3f",df[ln,args[i]])
                elseif precision == 4
                    nstr = @sprintf("%.4f",df[ln,args[i]])
                elseif precision == 5
                    nstr = @sprintf("%.5f",df[ln,args[i]])
                elseif precision == 6
                    nstr = @sprintf("%.6f",df[ln,args[i]])
                elseif precision == 7
                    nstr = @sprintf("%.7f",df[ln,args[i]])
                elseif precision == 8
                    nstr = @sprintf("%.8f",df[ln,args[i]])
                else
                    error("Can't display that much precision.")
                end
                print(prepend_spaces(nstr,lenvec[i]))
            end
            if i < nargs
                print("  ")
            end
        end
        print("\n")
    end
end


function desc(df::DataFrame;label_dict::Dict=Dict())

    varnames = names(df)
    maxval = maximum(map(x->haskey(label_dict["label"],string(x)) ? length(string(x)) : 0, varnames))
    maxval = maxval < 8 ? 8 : maxval
    maxtype = 7
    maxlab = maximum(map(x->haskey(label_dict["label"],string(x)) ? length(label_dict["label"][string(x)]) : 0, varnames))
    maxlab = maxlab < 11 ? 11 : maxlab
    maxformat = maximum(map(x->haskey(label_dict["label"],string(x)) ? length(label_dict["format"][string(x)]) : 0, varnames))
    maxformat = maxformat < 6 ? 6 : maxformat

    numvar = length(varnames)

    # variable index
    maxobs = length(string(numvar))
    print(prepend_spaces("Num",maxobs < 3 ? 3 : maxobs),"  ")

    # variable name
    print(append_spaces("Variable",maxval),"  ")

    # type (eltype)
    print(append_spaces("Type",maxtype),"  ")

    if isempty(label_dict) == false
        # value label
        print(append_spaces("Value Label",maxlab),"  ")

        # format
        print(append_spaces("Format",maxformat),"  ")

        # label
        print("Label")
    end

    print("\n")

    # dashes for a line by itself -- assume 20 characters for "label"
    if isempty(label_dict) == true
        println(repeat("-",maxobs+maxval+maxtype+4))
    else
        println(repeat("-",maxobs+maxval+maxtype+maxlab+maxformat+30))
    end

    i=1
    for v::Symbol in names(df)
        eltyp = eltype(df[v])

        if eltyp <: AbstractString
            eltyp = "String"
        end

        v_str = string(v)
        print(prepend_spaces(string(i),maxobs),"  ",append_spaces(v_str,maxval),"  ",append_spaces(string(eltyp),maxtype),"  ")

        if haskey(label_dict["value"],v_str)
            print(append_spaces(label_dict["value"][v_str],maxlab),"  ")
        else
            print(repeat(" ",maxlab+2))
        end

        if haskey(label_dict["format"],v_str)
            print(append_spaces(label_dict["format"][v_str],maxformat),"  ")
        else
            print(repeat(" ",maxformat+1))
        end

        if haskey(label_dict["variable"],v_str)
            print(label_dict["variable"][v_str])
        end

        print("\n")
        i += 1
    end
end

function getmaxlength(da::DataArray)
    maxlen = 0

    for i = 1:length(da)
        if isna(da[i])
            continue
        end
        maxlen = max(maxlen,length(da[i]))
    end

    return maxlen
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
                if args[1] <= maxlength
                    push!(dslist,v)
                end
            elseif length(args) == 2
                if args[1] <= maxlength <= args[2]
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


import DataFrames.dropna

"""
    dropna(df::DataFrame,args::Symbol...)

Create a data frame that does not contain NA values in `args` columns. When no `args`
are specified, it produces a data frame that does not not contain any NA values
(e.g., it is identical to `df = df[completecases(df),:]`).
"""
function dropna(df::DataFrame,args::Symbol...)
    if length(args) == 0
        return df[completecases(df),:]
    end

    ba = BitArray(trues(size(df,1)))
    for sym in args
        ba &= ~BitArray(df[sym].na)
    end
    return df[ba,:]
end

"""
    andcond(df::DataFrame,args::Expr...)

Produces a BitArray array in which rows that meet all the `args` conditions will have
a `false` value. It is intended to be used to drop cases from a DataFrame.
Each condition in `args` must be be an Expr: e.g, `:(condition)`.
"""
function andcond(df::DataFrame,args::Expr...)
    if length(args) == 0
        error("No condition was specified.")
    end

    ba = BitArray(trues(size(df,1)))
    for arg in args
        ba &= BitArray(eval(arg))
    end
    return ba
end

"""
    orcond(df::DataFrame,args::Expr...)

Produces a BitArray array in which rows that meet at least one `args` condition will have
a `false` value. It is intended to be used to drop cases from a DataFrame.
Each condition in `args` must be an Expr: e.g, `:(condition)`.
"""
function orcond(df::DataFrame,args::Expr...)
    if length(args) == 0
        error("No condition was specified.")
    end

    ba = BitArray(trues(size(df,1)))
    for arg in args
        ba |= BitArray(eval(arg))
    end
    return ba
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
#duplicates(df::DataFrame,arg::Symbol; cmd::Symbol = :report) = duplicates(df,Tuple(arg),cmd=cmd)

"""
    sampleselect(df::DataFrame,num::Union{Int64,Float64})

Creates a dataframe with a randomly selected sample from the input dataframe. `num`
specifies the sample size. If `num` is an integer, it indicates the number of rows to be selected.
If it is a float, it indicates the percentage of the input dataframe to be randomly selected.

##Example
To select 100 rows, use

```jldoctest
julia> df2 = sampleselect(df,100)
```

To select 20% of the original dataframe, use

```jldoctest
julia> df2 = sampleselect(df,20.)
```

or

```jldoctest
julia> df2 = sampleselect(df,.2)
```

"""
function sampleselect(df::DataFrame,num::Union{Int64,Float64})

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
