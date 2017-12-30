"""
    tab(df::DataFrame,vars::Symbol...; label_dict::Union{Void,Dict} = nothing, rmna::Bool = true, weights::AbstractVector = UnitWeights())

Produce n-way frequency table from a DataFrame or any type of arrays. Weights can be used to obtain
weighted frequencies. `tab` is mainly a wrapper for the excellent `FreqTables` package for all
but DataArrays with NA values. Use `rmna = false` to obtain frequencies that include NA values.
The returned table is in `NamedArrays`. Frequencies are in an n-dimensional array `na.array`
where `na` is the returned NamedArray. `na.dimnames` contain row and column values.
When NA values are included, the `dimnames` (see `NamedArrays` package) are returned in string
values rather than in the original data type.
"""
function tab(df::DataFrame,args::Symbol...; rmna = true, weights::AbstractVector = FreqTables.UnitWeights())

    # number of complete cases
    ba = completecases(df[collect(args)])
    cc = sum(ba)

    # remove NAs
    if rmna == true
        df = df[ba,collect(args)]
    end

    # find out if the args columns contain any NA values
    if rmna == true || cc == size(df,1)
        a = FreqTables.freqtable([df[ba,y] for y in args]...; weights = weights)
    else
        # there are NA values and so we cannot use freqtable
        # weights are not allowed, either
        a = tabna([df[y] for y in args]...)
    end

    setdimnames!(a, collect(args))

    return a
end
function tabna(x::AbstractVector...)

    ncols = length(x)
    nrows = length(x[1])
    for i = 2:ncols
        if nrows != length(x[i])
            error("Columns do not have the same rows.")
        end
    end

    # create output arrays
    vdims = Vector{Vector}(ncols)
    vnums = zeros(Int64,ncols)
    naidx = falses(ncols)
    for i = 1:ncols
        if sum(x[i].na) > 0 # there are NA values in this vector
            vdims[i] = sort(levels(dropna(x[i])))
            vnums[i] = size(vdims[i],1) + 1 # one for the NA
            naidx[i] = true
        else
            vdims[i] = sort(levels(x[i]))
            vnums[i] = size(vdims[i],1)
        end
    end

    # allocate memory for the output array
    a = zeros(Int64,vnums...)
    idxvec = zeros(Int64,ncols)

    # get frequencies with the NA values in the vnums... cell for each column
    for i = 1:nrows
        for j = 1:ncols
            @inbounds idxvec[j] = isna(x[j][i]) ? vnums[j] : searchsortedfirst(vdims[j],x[j][i])
        end
        @inbounds a[idxvec...] += 1
    end

    dimnames = Vector{Vector}(ncols)
    for i = 1:ncols
        if naidx[i]
            dimnames[i] = vcat(string.(vdims[i]),"NA")
        else
            dimnames[i] = vdims[i]
        end
    end

    return NamedArray(a, tuple(dimnames...), ntuple(i -> "Dim$i", ncols))
end
