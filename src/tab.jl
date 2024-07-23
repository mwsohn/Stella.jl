"""
    tab(df::DataFrame,vars::Symbol...; labels = nothing, maxrows = -1, maxcols = 20, decimals=4)
    tab(na::NamedArray; labels = nothing, maxrows = -1, maxcols = 20, decimals=4)

Produce an one-way or two-way frequency table from a DataFrame or a NamedArray obtained from
freqtable function. `tab` is mainly a wrapper for the excellent `FreqTables` package.
Use `skipmissing = true` to obtain frequencies that include `missing` values.
The returned table is a `NamedArray`. Frequencies are in an n-dimensional array `na.array`
where `na` is the returned NamedArray. 
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
function tab(indf,var::Union{Symbol,String}; decimals=4, skipmissing=true)
    if in(string(var),names(indf)) == false
        println("$var is not found in the input DataFrame.")
        return nothing
    end
    _tab1(freqtable(indf,var, skipmissing=skipmissing); decimals=decimals)
end
function tab(indf,var1::Union{Symbol,String},var2::Union{Symbol,String}; maxrows = -1, maxcols = 20, decimals=4, skipmissing=true)
    if in(string(var1),names(indf)) == false
        println("$var1 is not found in the input DataFrame.")
        return nothing
    end
    if in(string(var2), names(indf)) == false
        println("$var2 is not found in the input DataFrame.")
        return nothing
    end
    _tab2(freqtable(indf, var1, var2, skipmissing=skipmissing); maxrows=maxrows, maxcols=maxcols, decimals=decimals)
end
function tab(indf,var1::Union{Symbol,String},var2::Union{Symbol,String},var3::Union{Symbol,String};
    maxrows=-1, maxcols=20, decimals=4, skipmissing=true)
    for v in (var1,var2,var3)
        if in(string(v), names(indf))
            println("$v is not found in the input DataFrame.")
            return nothing
        end
    end
    _tab3(freqtable(indf, var1, var2, var3, skipmissing=skipmissing); maxrows=maxrows, maxcols=maxcols, decimals=decimals)
end

function _tab1(na::NamedArray; decimals = 4)
 
    # do not output rows with zeros
    z = findall(x -> x != 0, na.array)
    arry = na.array[z]

    # value labels and "Total"
    rownames = vcat(names(na)[1][z],"Total")

    # counts - the last row has the total
    counts = vcat(arry,sum(na,dims=1))

    # percents
    percents = 100 .* counts ./ counts[end]

    # cumulative percents
    cumpct = 100 .* vcat(cumsum(arry,dims=1),counts[end]) ./ counts[end]

    ar = hcat(rownames,counts, percents, cumpct)

    PrettyTables.pretty_table(ar, 
        header=[na.dimnames[1],"Counts"," Percent","Cum Pct"],
        formatters = ft_round(decimals,[3,4]),
        crop = :none,
        hlines=[0,1,length(rownames),length(rownames)+1],
        vlines=[1])
end

function _tab2(na::NamedArray; maxrows = -1, maxcols = 20, decimals=4)
    
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
    (statistic, dof, pval) = Stella.chi2(testarray)
    println("Pearson chi-square = ", @sprintf("%.4f",statistic), " (", dof, "), p ", 
        pval < 0.0001 ? "< 0.0001" : string("= ",round(pval,sigdigits = 6)))

    if size(testarray) == (2, 2) # 2x2 array
        println("Fisher's exact test = ", @sprintf("%.4f",
            pvalue(HypothesisTests.FisherExactTest(rec(testarray')...))))
    end
end

function _tab3(na::NamedArray; maxrows = -1, maxcols = 20, decimals=4)

    # stratify the var3 (na.dimnames[3])
    n3 = size(na,3)
    vals = na.dicts[3].keys

    for i in 1:n3
        println("\n\n",na.dimnames[3], " = ", vals[i] ,"\n")

        _tab2(na[:,:,i]; maxrows = maxrows, maxcols = maxcols, decimals=decimals)
    end
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
