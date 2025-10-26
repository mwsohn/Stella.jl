"""
    tab(df::DataFrame,vars::Symbol...; maxrows = -1, maxcols = 20, sort = true, skipmissing=true)
    tab(na::NamedArray; maxrows = -1, maxcols = 20)

Produce an one-way or two-way frequency table from a DataFrame or a NamedArray obtained from
freqtable function. `tab` is mainly a wrapper for the excellent `FreqTables` package.

Use `skipmissing = false` to obtain frequencies that include `missing` values.

For an one-way frequency table, the table can be
ordered by the frequency by specifying `sort = true` as an option.

For a two-way table, summary values of a continuous variable can be requested by specifying the variable name
as an option `summrize = :var`.
"""
function tab(na::NamedArray; skipmissing=true)
    
    len = length(na.dimnames)
    if len == 1
        _tab1(na; skipmissing=skipmissing)
    elseif len == 2
        _tab2(na; skipmissing=skipmissing)
    elseif len == 3
        _tab3(na; skipmissing=skipmissing)
    else
        error("Crosstabs for more than 3 variables are not currently supported.")
    end
end
function tab(indf,var::Union{Symbol,String}; skipmissing=true, sort=false, summarize=nothing)
    if in(string(var),names(indf)) == false
        throw(ArgumentError("$var is not found in the input DataFrame."))
        return nothing
    end
    if summarize != nothing
        _tab1summarize(indf,var,summarize)
    end
    _tab1(freqtable(indf,var, skipmissing=skipmissing); sort=sort)
end
function tab(ivar::AbstractVector; skipmissing=true, sort=false)
    _tab1(freqtable(ivar, skipmissing=skipmissing); sort=sort)
end
function tab(indf,var1::Union{Symbol,String},var2::Union{Symbol,String}; maxrows = -1, maxcols = 20, skipmissing=true, summarize = nothing)
    if summarize != nothing && isa(indf[:,summarize], CategoricalArray)
        throw(ArgumentError("$summarize cannot be a CategoricalArray."))
        return nothing
    end
    if in(string(var1),names(indf)) == false
        println("$var1 is not found in the input DataFrame.")
        return nothing
    end
    if in(string(var2), names(indf)) == false
        println("$var2 is not found in the input DataFrame.")
        return nothing
    end
    if summarize == nothing
        return _tab2(freqtable(indf, var1, var2, skipmissing=skipmissing); maxrows=maxrows, maxcols=maxcols)
    end

    _tab2summarize(indf, var1, var2, summarize; maxrows=-1, maxcols=20)
end
function tab(indf,var1::Union{Symbol,String},var2::Union{Symbol,String},var3::Union{Symbol,String};
    maxrows=-1, maxcols=20, skipmissing=true, summarize=nothing)
    for v in (var1,var2,var3)
        if in(string(v), names(indf)) == false
            println("$v is not found in the input DataFrame.")
            return nothing
        end
    end
    
    if summarize == nothing
        na = freqtable(indf, var1, var2, var3, skipmissing=skipmissing)
        # stratify the var3 (na.dimnames[3])
        n3 = size(na, 3)
        vals = na.dicts[3].keys

        for i in 1:n3
            println("\n\n", na.dimnames[3], " = ", vals[i], "\n")

            _tab2(na[:, :, i]; maxrows=maxrows, maxcols=maxcols)
        end
    else
        n3 = sort(unique(indf[!,var3]))
        for v in n3
            println("\n\n", var3, " = ", v, "\n")

            _tab2summarize(filter(x -> x[var3] == v, indf), var1, var2, summarize; maxrows=-1, maxcols=20)
        end
    end
end

function tabi(a::AbstractArray)
    if length(size(a)) == 2 && all(x -> x >= 2,a)
        Stella._tab2(NamedArray(a))
    end
    error("Input array must be 2x2 and have at least two levels on each dimension.")
end

function _tab1(na::NamedArray; sort = false)
 
    # do not output rows with zeros
    z = findall(x -> x != 0, na.array)
    arry = na.array[z]

    if sort
        s = sortperm(arry,rev=true)
        arry = arry[s]
    end

    # value labels and "Total"
    if sort
        rownames = vcat(names(na)[1][z][s],"Total")
    else
        rownames = vcat(names(na)[1][z],"Total")
    end

    # counts - the last row has the total
    counts = vcat(arry,sum(na,dims=1))

    # percents
    percents = 100 .* counts ./ counts[end]

    # cumulative percents
    cumpct = 100 .* vcat(cumsum(arry,dims=1),counts[end]) ./ counts[end]

    ar = hcat(rownames,counts, percents, cumpct)

    PrettyTables.pretty_table(ar, 
        header=[na.dimnames[1],"Counts"," Percent","Cum Pct"],
        formatters = ft_round(3,[3,4]),
        crop = :none,
        hlines=[0,1,length(rownames),length(rownames)+1],
        vlines=[1])
end

function _tab2(na::NamedArray; maxrows = -1, maxcols = 20)
  
    # counts
    counts = na.array
    counts = vcat(counts,sum(counts,dims=1)) # column sum
    counts = hcat(counts,sum(counts,dims=2)) # row sum
    rz = findall(x -> x != 0, counts[:, end]) # find all columns with non-zero totals
    cz = findall(x -> x != 0, counts[end, :]) # find all rows with non-zero totals
    counts = counts[rz,cz]
    (nrow, ncol) = size(counts)

    # value labels and "Total"
    rownames = vcat(names(na)[1], "Total")[rz]

    # colunm names
    colnames = vcat(names(na)[2], "Total")[cz]

    # row and column percentages
    rowpct = 100 .* counts ./ counts[:,end]
    colpct = (100 .* counts' ./ counts[end,:])'

    # interleave them 
    d = reshape(Any[counts rowpct colpct]'[:],(ncol,(nrow)*3))'

    # add two blank cells
    rownames2 = vcat([ [x, " ", " "] for x in rownames ]...)


    pretty_table(d,
        row_labels = rownames2,
        row_label_column_title=string(na.dimnames[1], " / ", na.dimnames[2]),
        header=colnames,
        crop = :none, 
        max_num_of_rows = maxrows,
        max_num_of_columns = maxcols,
        hlines=vcat([0, 1], [x * 3 + 1 for x in 1:(nrow+1)]),
        vlines = [1])

    testarray = na.array[rz[1:end-1],cz[1:end-1]]
    if size(testarray,1) > 1 && size(testarray,2) > 1
        (statistic, dof, pval) = Stella.chi2(testarray)
        println("Pearson chi-square = ", @sprintf("%.4f",statistic), " (", dof, "), p ", 
            pval < 0.0001 ? "< 0.0001" : string("= ",round(pval,sigdigits = 6)))
    end

    if size(testarray) == (2, 2) && all(x -> x > 0, testarray) # 2x2 array
        println("Fisher's exact test = ", @sprintf("%.4f",
            pvalue(HypothesisTests.FisherExactTest((testarray')...))))
    end
end

function _tab1summarize(indf,var,sumvar)
    ba = completecases(indf[!, [var, sumvar]])
    odf = combine(groupby(indf[ba,[var,sumvar]],var), nrow => :n, sumvar => mean => :mean, sumvar => std => :sd)
    tdf = DataFrame(n = nrow(indf[ba,:]), mean = mean(indf[ba,sumvar]), sd = std(indf[ba,sumvar]))
    odf = vcat(odf,tdf)
    pretty_table(odf,
        row_labels=odf[!,sumvar],
        row_label_column_title=string(sumvar),
        header=["N","Mean","StDev"],
        crop=:none,
        max_num_of_rows=maxrows,
        max_num_of_columns=maxcols,
        hlines=vcat([0, 1], nrows(odf)+1),
        vlines=[1])
end

function _tab2summarize(indf, var1, var2, sumvar; maxrows=-1, maxcols=20)
    ba = completecases(indf[!,[var1,var2,sumvar]])
    na = freqtable(indf[ba,:], var1, var2)

    # margin stats
    gdf = groupby(indf[ba, :], var1)
    var1df = combine(gdf, nrow => :n, sumvar => mean => :mean, sumvar => std => :sd)
    nrows = size(var1df,1)

    gdf = groupby(indf[ba, :], var2)
    var2df = combine(gdf, nrow => :n, sumvar => mean => :mean, sumvar => std => :sd)
    ncols = size(var2df,1)

    # cell stats
    gdf = groupby(indf[ba,:],[var1,var2])
    outdf = combine(gdf, nrow => :n, sumvar => mean => :mean, sumvar => std => :sd)

    # value labels and "Total"
    rownames = vcat(collect(skipmissing(names(na)[1])), "Total")
    rownames2 = vcat([[x, " ", " "] for x in rownames]...)

    # colunm names
    colnames = vcat(collect(skipmissing(names(na)[2])), "Total")

    # combine stats
    d = Any[outdf.mean outdf.sd outdf.n]'
    for i = 1:nrows
        idx = (i - 1) * ncols + 1
        if i == 1
            e = d[1:3, 1:ncols]
        else
            e = vcat(e, d[1:3, idx:idx+ncols-1])
        end
    end

    # row margins
    e = hcat(e, Any[var1df.mean var1df.sd var1df.n]'[:])

    # column margins
    cm = Any[var2df.mean var2df.sd var2df.n]'[:]

    # grand total
    push!(cm, mean(indf[ba,sumvar]))
    push!(cm, std(indf[ba, sumvar]))
    push!(cm, size(indf[ba, sumvar],1))

    cm = reshape(cm,(3,ncols + 1))

    # combine cell summary stats with column margin stats
    e = vcat(e,cm)

    # output
    pretty_table(e,
        row_labels=rownames2,
        row_label_column_title=string(na.dimnames[1], " / ", na.dimnames[2]),
        header=colnames,
        crop=:none,
        max_num_of_rows=maxrows,
        max_num_of_columns=maxcols,
        hlines=vcat([0, 1], [x * 3 + 1 for x in 1:(nrows + 1)]),
        vlines=[1])

end

"""
    chi2(m::AbstractMatrix{Integer})

Returns chi squared test statistic and p-value for the test of independence.
"""
function chi2(m::AbstractMatrix{T}) where {T <: Integer}
    (nrow, ncol) = size(m)
    if nrow <= 1 || ncol <= 1
        error("at least a 2x2 table is expected")
    end
    rowsum = sum(m, dims=2)
    colsum = sum(m, dims=1)
    dof = (nrow - 1)*(ncol - 1)
    e = rowsum * colsum ./ sum(m)
    statistic = sum((m .- e).^2 ./ e)
    pvalue = Distributions.ccdf(Distributions.Chisq(dof),statistic)
    return (statistic, dof, pvalue)
end
