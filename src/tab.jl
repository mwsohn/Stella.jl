"""
    tab(df::DataFrame,vars::Symbol...; labels = nothing, maxrows = -1, maxcols = 20, decimals=4)
    tab(na::NamedArray; labels = nothing, maxrows = -1, maxcols = 20, decimals=4)

Produce an one-way or two-way frequency table from a DataFrame or a NamedArray obtained from
freqtable function. `tab` is mainly a wrapper for the excellent `FreqTables` package.
Use `skipmissing = true` to obtain frequencies that include `missing` values.
The returned table is a `NamedArray`. Frequencies are in an n-dimensional array `na.array`
where `na` is the returned NamedArray. 
"""
function tab(na::NamedArray)
    
    len = length(na.dimnames)
    if len == 1
        _tab1(na)
    elseif len == 2
        _tab2(na) 
    elseif len == 3
        _tab3(na)
    else
        error("Crosstabs for more than 3 variables are not currently supported.")
    end
end
function tab(indf,var::Union{Symbol,String}; decimals=3, labels=nothing)
    _tab1(freqtable(indf,var); decimals=decimals, labels=labels)
end
function tab(indf,var1::Union{Symbol,String},var2::Union{Symbol,String}; 
    maxrows = -1, maxcols = 20, labels=nothing)
    _tab2(freqtable(indf,var1,var2); maxrows=maxrows, maxcols = maxcols, labels = labels)
end
function tab(indf,var1::Union{Symbol,String},var2::Union{Symbol,String},var3::Union{Symbol,String};
    maxrows = -1, maxcols = 20, labels=nothing)
    _tab3(freqtable(df,var1,var2,var3); maxrows=maxrows, maxcols=maxcols, labels=labels)
end


function _tab1(na::NamedArray; decimals = 4, labels=labels)
 
    # rownames
    rownames = names(na)[1]

    # Total row label
    rownames = vcat(rownames,"Total")

    # counts
    counts = vcat(na.array,sum(na,dims=1))

    # percents
    percents = 100 .* counts ./ counts[end]

    # cumulative percents
    cumpct = 100 .* vcat(cumsum(na.array,dims=1),counts[end]) ./ counts[end]

    ar = hcat(rownames,counts, percents, cumpct)

    PrettyTables.pretty_table(ar, 
        header=[na.dimnames[1],"Counts","Percent","Cum. Percent"],
        formatters = ft_round(decimals,[3,4]),
        hlines=[0,1,length(rownames),length(rownames)+1])
end

function _tab2(na::NamedArray; maxrows = -1, maxcols = 20, labels=nothing)
    
    # rownames
    rownames = names(na)[1]
    rownames = vcat(rownames,"Total")

    # colunm names
    colnames = names(na)[2]
    colnames = vcat(colnames,"Total") 

    # counts
    counts = na.array
    (nrow,ncol) = size(counts)
    counts = vcat(counts,sum(counts,dims=1)) # column sum
    counts = hcat(counts,sum(counts,dims=2)) # row sum

    # row and column percentages
    rowpct = 100 .* counts ./ counts[:,ncol+1]
    colpct = (100 .* counts' ./ counts[nrow+1,:])'

    # interleave them 
    d = reshape(Any[counts rowpct colpct]'[:],(ncol+1,(nrow+1)*3))'

    # add two blank cells
    rownames2 = vcat([ [x, " ", " "] for x in rownames ]...)

    # labels
    # if labels != nothing && varlab(labels,na.dimnames[1]) != nothing
        

    # hline location
    hlines = vcat([0,1], [ x*3 + 1 for x in 1:(nrow+1) ])


    pretty_table(hcat(rownames2,d),
        header=vcat( string(na.dimnames[1]," / ", na.dimnames[2]), colnames),
        # formatters = ft_round(decimals),
        max_num_of_rows = maxrows,
        max_num_of_columns = maxcols,
        hlines = hlines)

    (statistic, pval) = Stella.chi2(na.array)
    println("Pearson chi-square = ", @sprintf("%.4%",statistic), " (", (ncol-1)*(nrow-1), "), 
        p ", pval < 0.0001 ? "< 0.0001" : string("= ",round(pval,sigdigits = 6)))
end

function _tab3(na::NamedArray; maxrows = -1, maxcols = 20, labels=nothing)

    # stratify the var3 (na.dimnames[3])
    n3 = size(na,3)
    vals = na.dicts[3].keys

    for i in 1:n3
        println("\n\n",na.dimnames[3], " = ", vals[i] ,"\n")

        _tab2(na[:,:,i]; maxrows = maxrows, maxcols = maxcols, labels=labels)
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
    df = (nrow - 1)*(ncol - 1)
    e = rowsum * colsum ./ sum(m)
    statistic = sum((m .- e).^2 ./ e)
    pvalue = HypothesisTests.ccdf(Distributions.Chisq(df),statistic)
    return (statistic, pvalue)
end
