"""
    tab(df::DataFrame,vars::Symbol...; labels::Union{Void,Label} = nothing, skipmissing::Bool = true, weights::AbstractWeights = fweights(ones(Int8,length(vars[1]))))
    tab(a::AbstractArray...; skipmissing::Bool = true, weights::AbstractWeights = fweights(ones(Int8,length(vars[1]))))

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
    @inbounds @simd for i in 1:alen
        index = Base.ht_keyindex(d, A[i])

        if index > 0
            @inbounds d.vals[index] += weights[i]
        else
            @inbounds d[A[i]] = weights[i]
        end
    end

    # construct labels
    alev = sort(collect(keys(d)))

    cnt = Vector{Int64}(length(alev))
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

    if nvec == 1
        na = tab(df[vars[1]], skipmissing = skipmissing, weights = weights)
        setdimnames!(na,(vars[1],"Stat"))
        return na
    end
    na = tab(df[vars[1]],df[vars[2]], skipmissing = skipmissing, weights = weights)
    # if labels != nothing # && Labels.defined(labels,v)
    #     # value labels
    #     labval = Vector{String}(length(vec))
    #     for (i,val) in enumerate(lev)
    #         labval[i] = vallab(labels,v,val)
    #     end
    #
    #     # allowmissing
    #     if allowmissing == true
    #         push!(labval,"missing")
    #     end
    #
    #     # variable label
    #     labvar = varlab(labels,v)
    #
    #     # return a NamedArray
    #     setnames!(na,labval,1)
    #     setnames!(na,["Frequency"],2)
    #     setdimnames!(na,[labvar,"count"])
    #     return na
    # end

    setdimnames!(na,(vars[1],vars[2]))
    return na
end

struct XsqResult
    chisq::Float64
    dof::Int
    p::Float64
end
function chi2test(a::Array{Float,2})

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

function allpercent(tab::NamedArray)

    totrow = sum(tab,2)
    totcol = sum(tab,1)
    total  = sum(tab)

    return (NamedArray(100 * tab ./ totrow,tuple(names(tab)...),tuple(dimnames(tab)...)),
        NamedArray(100 * tab ./ totcol,tuple(names(tab)...),tuple(dimnames(tab)...)),
        NamedArray(100 * tab ./ total,tuple(names(tab)...),tuple(dimnames(tab)...)))
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
