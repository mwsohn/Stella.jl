"""
    tab(df::DataFrame,vars::Symbol...; label_dict::Union{Nothing,Dict} = nothing, rmna::Bool = true, weights::AbstractVector = UnitWeights())

Produce an one-way or two-way frequency table from a DataFrame or any type of arrays.
Any weights defined under AbstractWeights (e.g., fweights, aweights, or pweights)
can be used to obtain weighted frequencies.

`tab` is mainly a wrapper for the excellent `FreqTables` package.
Use `skipmissing = true` to obtain frequencies that include `missing` values.
The returned table is a `NamedArray`. Frequencies are in an n-dimensional array `na.array`
where `na` is the returned NamedArray. `na.dimnames` contain row and column values.
"""
function tab(A::AbstractArray,B::AbstractArray; skipmissing::Bool=false,weights::AbstractWeights=fweights(ones(Int8,length(A))))

    # check lengths of arrays
    alen = length(A)
    blen = length(B)
    wlen = length(weights)

    if alen != blen
        error("Lengths of two vectors are not equal")
    end

    if alen != wlen
        error("Input vectors and `weights` vector are not of the same lengths")
    end

    # from FreqTables.jl/src/freqtable.jl
    d = Dict{Tuple{eltype(A),eltype(B)},Int64}()
    @inbounds @simd for i in 1:alen

        el = (A[i],B[i])

        index = Base.ht_keyindex(d, el)

        if index > 0
            @inbounds d.vals[index] += weights[i]
        else
            @inbounds d[el] = weights[i]
        end
    end

    # construct labels
    alev = Set{eltype(A)}()
    blev = Set{eltype(B)}()
    for (k1,k2) in keys(d)
        push!(alev,k1)
        push!(blev,k2)
    end

    alev = sort!(collect(alev))
    blev = sort!(collect(blev))
    if skipmissing
        alev = collect(skipmissing(alev))
        blev = collect(skipmissing(blev))
    end

    cnt = Array{Int64,2}(length(alev),length(blev))
    @inbounds @simd for i=1:length(alev)
        @inbounds @simd for j=1:length(blev)
            el = (alev[i],blev[j])
            if haskey(d,el)
                cnt[i,j] = d[el]
            else
                cnt[i,j] = 0
            end
        end
    end

    return NamedArray(cnt,(alev,blev),("A","B"))
end
function tab(A::AbstractArray; skipmissing::Bool=true,weights::AbstractWeights=fweights(ones(Int8,length(A))))

    alen = length(A)
    wlen = length(weights)

    if alen != wlen
        error("The `weights` vector does not have the same length as the input vector `A`")
    end

    # from FreqTables.jl/src/freqtable.jl
    d = Dict{eltype(A),Int64}()
    @simd for i in 1:alen
        index = Base.ht_keyindex(d, A[i])

        if index > 0
            @inbounds d.vals[index] += weights[i]
        else
            @inbounds d[A[i]] = weights[i]
        end
    end

    # construct labels
    alev = sort(collect(keys(d)))

    cnt = Vector{Int64}(undef,length(alev))
    @inbounds @simd for i=1:length(alev)
        cnt[i] = d[alev[i]]
    end

    if skipmissing == true && alev[end] == missing
        return NamedArray(reshape(cnt[1:end-1],length(alev)-1,1),(alev[1:end-1],["Frequency"]),("A","Stat"))
    end

    return NamedArray(reshape(cnt,length(alev),1),(alev,["Frequency"]),("A","Stat"))
end
function tab(df::DataFrame,vars::Symbol...;  skipmissing::Bool=false,weights::AbstractWeights=fweights(ones(Int8,size(df,1))))
    nvec = length(vars)
    if nvec == 0
        error("At least one variable is required")
    elseif nvec > 2
        error("More than two variables are not supported")
    end

    # oneway table
    if nvec == 1
        na = tab(df[vars[1]], skipmissing = skipmissing, weights = weights)
        setdimnames!(na,(vars[1],"Stat"))
        return na
    end

    # two-way table
    na = tab(df[vars[1]], df[vars[2]], skipmissing = skipmissing, weights = weights)

    setdimnames!(na,(vars[1],vars[2]))
    return na
end

struct XsqResult
    chisq::Float64
    dof::Int
    p::Float64
end
function chi2test(a::Array{Float64,2})

    if ndims(t) != 2
        error("Only two dimensional arrays are supported")
    end

    rowsum = Array(sum(t,2))
    colsum = Array(sum(t,1))
    total = sum(t)

    ncol = length(colsum)
    nrow = length(rowsum)

    chisq = 0.
    for i = 1:nrow
        for j = 1:ncol
            expected = rowsum[i]*colsum[j]/total
            chisq += ((t[i,j] - expected)^2)/expected
        end
    end

    # degress of freedom
    df = (ncol-1)*(nrow-1)

    # return a tuple of chisq, df, p-value
    return XsqResult(chisq,df,Distributions.ccdf(Distributions.Chisq(df),chisq))
end

function rowpercent(tab::NamedArray)

    totrow = sum(tab,2)

    array = 100 * tab.array ./ totrow

    return NamedArray(array,tuple(names(tab)...),tuple(dimnames(tab)...))

end

function colpercent(tab::NamedArray)

    totcol = sum(tab,1)

    array = 100 * tab.array ./ tocol

    return NamedArray(array,tuple(names(tab)...),tuple(dimnames(tab)...))
end

function cellpercent(tab::NamedArray)

    total = sum(tab.array)

    array = 100 * tab.array ./ total

    return NamedArray(array,tuple(names(tab)...),tuple(dimnames(tab)...))
end


function tabprint(na::NamedArray; precision=2, chisq=true, row=true, col=true, cell=false, all=false, pagewidth=78)

    if ndims(na) == 2 && length(names(na,2)) == 1 && names(na,2)[1] == "Frequency"
        tabprint1(na, all = all, precision = precision)
        exit()
    elseif ndims(na) > 2
        error("Only up to two dimensional arrays are currently supported")
    end

    print(dimnames(na,1), " \\ ", dimnames(na,2),"\n")

    # all is true when row, col, or cell is true
    if all
        row = col = cell = true
    end

    # row names
    rownames = string.(names(na,1))
    maxrowname = maximum(vcat(5,length.(rownames))) # minimum width is 5

    # column names and their widths
    colnames = string.(names(na,2))
    maxcolname = maximum(vcat(3,length.(colnames))) # minimum width is 3

    # width of data columns - the same as the greater of the length of the grand total
    # and 4 + precision
    tot = sum(na.array) # grand total
    colwidth = max(length(digits(Int(floor(tot)))),4+precision)

    # number of columns
    ncols = length(colnames)

    # number of rows
    nrows = length(rownames)

    # floating point numbers with three digits after decimal point
    if eltype(na.array) <: AbstractFloat
        colwidth += 3
    end

    # determine column widths
    colwidth = max(maxcolname,colwidth)

    #---------------------------------------------------
    # column totals
    colsum = sum(na.array,1)

    # row totals
    rowsum = sum(na.array,2)

    #---------------------------------------------------
    # determine how many columns can be printed within
    # the pagewidth
    maxcol = Int(floor((pagewidth - maxrowname) / colwidth))
    niter = Int(ceil(ncols / maxcol))
    nlastrow = ncols % maxcol

    for i in 1:niter

        #---------------------------------------------------
        # print header
        print(repeat(" ",maxrowname)," |")

        for j = 1:maxcol
            c = Int(maxcol*(i-1))+j
            if c > ncols
                break
            end
            print(" ",lpad(string(colnames[c]),colwidth))
        end

        print(" | ",lpad("Total",colwidth),"\n")

        print(repeat("-",maxrowname),"-+-",repeat("-",(colwidth+1)*(i == niter ? nlastrow : maxcol)),"+-",repeat("-",colwidth),"\n")

        #----------------------------------------------------
        # print values
        for r = 1:nrows

            # row name
            print(rpad(string(rownames[r]),maxrowname)," |")

            for j = 1:maxcol
                c = Int(maxcol*(i-1))+j
                if c > ncols
                    break
                end
                val = na.array[r,c]
                print(" ",lpad(val,colwidth))
            end

            # row total
            print(" |")

            val = rowsum[r]
            print(" ",lpad(val,colwidth),"\n")

            # row percentages
            if row
                print(repeat(" ",maxrowname)," |")
                for j = 1:maxcol
                    c = Int(maxcol*(i-1)) + j
                    if c > ncols
                        break
                    end
                    val = strval(100 * na.array[r,c] / rowsum[r],precision)
                    print(" ",lpad(val,colwidth," "))
                end

                # row percentage
                print(" |")

                val = strval(100.0,precision)
                print(" ",lpad(val,colwidth),"\n")
            end

            # column percentages
            if col
                print(repeat(" ",maxrowname)," |")
                for j = 1:maxcol
                    c = Int(maxcol*(i-1)) + j
                    if c > ncols
                        break
                    end
                    val = strval(100 * na.array[r,c] / colsum[c],precision)
                    print(" ",lpad(val,colwidth))
                end

                # column percent
                print(" |")

                val = strval(100 * rowsum[r] / tot,precision)
                print(" ",lpad(val,colwidth),"\n")
            end

            # cell percentages
            if cell
                print(repeat(" ",maxrowname)," |")
                for j = 1:maxcol
                    c = Int(maxcol*(i-1)) + j
                    if c > ncols
                        break
                    end
                    val = strval(100 * na.array[r,c] / tot,precision)
                    print(" ",lpad(val,colwidth))
                end

                # column percent
                print(" |")

                val = strval(100 * rowsum[r] / tot,precision)
                print(" ",lpad(val,colwidth),"\n")
            end

            if row || col || cell
                print(repeat("-",maxrowname),"-+-",repeat("-",(colwidth+1)*(i == niter ? nlastrow : maxcol)),"+-",repeat("-",colwidth),"\n")
            end

        end

        #----------------------------------------------------
        # Total
        if !(row || col || cell)
            print(repeat("-",maxrowname+1),"+",repeat("-",(colwidth+1)*(i == niter ? nlastrow : maxcol)),"-+-",repeat("-",colwidth),"\n")
        end
        print(rpad("Total",maxrowname," ")," |")

        for j = 1:maxcol
            c = Int(maxcol*(i-1)) + j
            if c > ncols
                break
            end
            val=colsum[c]
            print(" ",lpad(val,colwidth," "))
        end

        # Grand total
        val = strval(tot)
        print(" | ",lpad(val,colwidth," "),"\n")

        #----------------------------------------------------
        # row percentages
        if row
            print(repeat(" ",maxrowname)," |")
            for j = 1:maxcol
                c = Int(maxcol*(i-1)) + j
                if c > ncols
                    break
                end
                val = strval(100 * colsum[c] / tot,precision)
                print(" ",lpad(val,colwidth," "))
            end

            # column percent for the total column
            val = strval(100.0,precision)
            print(" | ",lpad(val,colwidth," "),"\n")
        end

        # column percentages
        if col
            print(repeat(" ",maxrowname)," |")
            for j = 1:maxcol
                c = Int(maxcol*(i-1)) + j
                if c > ncols
                    break
                end
                val = strval(100.0,precision)
                print(" ",lpad(val,colwidth," "))
            end

            # column percent for the total column
            val = strval(100.0,precision)
            print(" | ",lpad(val,colwidth," "),"\n")

        end

        # cell percentages
        if cell
            print(repeat(" ",maxrowname)," |")
            for j = 1:maxcol
                c = Int(maxcol*(i-1)) + j
                if c > ncols
                    break
                end
                val = strval(100.0 * colsum[c] / tot,precision)
                print(" ",lpad(val,colwidth," "))
            end

            # column percent
            val = strval(100.0,precision)
            print(" | ",lpad(val,colwidth," "),"\n")
        end

        # separator
        print("\n\n")
    end

    if chisq
        # chisq test
        dof = (nrows-1)*(ncols-1)
        chisqval = 0.
        for i = 1:nrows
            for j = 1:ncols
                expected = rowsum[i]*colsum[j]/tot
                chisqval += ((na.array[i,j] - expected)^2)/expected
            end
        end

        # degress of freedom
        dof = (ncols-1)*(nrows-1)

        # return a tuple of chisq, df, p-value
        pval = Distributions.ccdf(Distributions.Chisq(dof),chisqval)

        # chisquare output
        println("Pearson χ² (",dof,") = ",@sprintf("%.5f",chisqval)," Pr = ",@sprintf("%.5f",pval))
    end
end

function tabprint1(na::NamedArray; skipmissing = false, precision=2, labels::Union{Label,Nothing} = nothing)

    if ndims(na) == 2 && names(na,2)[1] != "Frequency"
        error("Cannot print 2 or higher dimensions")
    end

    # variable name
    dname = dimnames(na,1)

    # value names
    if labels == nothing
        rownames = string.(names(na,1))
    else
        rownames = [ ismissing(x) ? "missing" : vallab(labels,dname,x) for x in names(na,1) ]
    end

    # maximum width for row labels
    if skipmissing == true && rownames[end] == "missing"
        maxrowname = maximum(vcat(5,length(string(dname)),length.(rownames[1:end-1]))) # minimum width is 5
    else
        maxrowname = maximum(vcat(5,length(string(dname)),length.(rownames))) # minimum width is 5
    end

    # number of rows
    nrows = length(rownames)

    # column names
    colnames = ["Frequency","Percent","Cum Percent"]

    # grand total
    if skipmissing == true && rownames[end] == "missing"
        tot = sum(na.array[1:end-1])
    else
        tot = sum(na)
    end

    # maximum column width
    maxcolwidth = maximum(vcat(length(string(tot)),length.(colnames)))

    # header
    print("\n")
    print(lpad(dname,maxrowname)," |")
    for i = 1:length(colnames)
        print(lpad(colnames[i],maxcolwidth+1))
    end
    print("\n")

    # separator
    print(repeat("-",maxrowname),"-+")
    for i = 1:length(colnames)
        print(repeat("-",maxcolwidth + 1))
    end
    print("\n")

    # rows
    cumpct = 0.
    for i = 1:nrows

        if skipmissing == true && rownames[i] == "missing"
            continue
        end

        print(lpad(rownames[i],maxrowname)," |")

        # frequency
        print(lpad(na.array[i],maxcolwidth + 1))

        # percent
        pct = na.array[i] / tot
        print(lpad(strval(100*pct,precision),maxcolwidth + 1))

        # cumulative percent
        cumpct += pct
        println(lpad(strval(100*cumpct,precision),maxcolwidth + 1))
    end

    # separator
    print(repeat("-",maxrowname),"-+")
    for i = 1:length(colnames)
        print(repeat("-",maxcolwidth + 1))
    end
    print("\n")

    # total
    print(lpad("Total",maxrowname)," |")

    # total frequency
    print(lpad(tot,maxcolwidth + 1))

    # percent
    print(lpad(strval(100.0,precision),maxcolwidth + 1))

    # cumulative percent
    println(lpad(strval(100.0,precision),maxcolwidth + 1))

    print("\n")
end




# function tab(x::AbstractVector...; allowmissing=false,weights::AbstractWeights = fweights(ones(Int8,length(x[1]))))
#
#     ncols = length(x)
#     nrows = length(x[1])
#     for i = 2:ncols
#         if nrows != length(x[i])
#             error("Columns do not have the same rows.")
#         end
#     end
#
#     vectypes = Missings.T.(eltype.(x))
#     println(vectypes)
#
#     # # create output arrays
#     # vdims = Vector{Vector}(ncols)
#     # vnums = zeros(Int64,ncols)
#     # naidx = falses(ncols)
#     # for i = 1:ncols
#     #     if sum(x[i].na) > 0 # there are NA values in this vector
#     #         vdims[i] = sort(levels(dropna(x[i])))
#     #         vnums[i] = size(vdims[i],1) + 1 # one for the NA
#     #         naidx[i] = true
#     #     else
#     #         vdims[i] = sort(levels(x[i]))
#     #         vnums[i] = size(vdims[i],1)
#     #     end
#     # end
#     # from FreqTables.jl/src/freqtable.jl
#     d = Dict{Tuple{Missings.T(eltype(a)),Missings.T(eltype(b))},Int64}()
#     @inbounds @simd for i in 1:alen
#         el = (a[i],b[i])
#         index = Base.ht_keyindex(d, el)
#
#         if index > 0
#             @inbounds d.vals[index] += weights[i]
#         else
#             @inbounds d[el] = weights[i]
#         end
#     end
#
#     # output array
#     cnt = Array{Int64,2}(length(alev),length(blev))
#     @inbounds @simd for i=1:length(alev)
#         @inbounds @simd for j=1:length(blev)
#             cnt[i,j] = d[alev[i],blev[j]]
#         end
#     end
#
#
#     return NamedArray(a, tuple(dimnames...), ntuple(i -> "Dim$i", ncols))
# end
