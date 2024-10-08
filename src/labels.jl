@reexport import TableMetadataTools

function labels!(df::AbstractDataFrame,dd::Dict)
    for key in keys(dd)
        label!(df,key,dd[key])
    end
end

function values!(df::AbstractDataFrame,v::Union{Symbol,String},val::Union{Dict,NamedTuple}; ordered = true)

    if isa(df[:,v],CategoricalArray) 
        uncategorical!(df,v)
    end

    if nonmissingtype(eltype(df[:,v])) <: Integer
        vv = categorical(df[!,v], ordered = true)
        df[!,v] = recode(vv,val...)
        levels!(df[!,v], last.(sort(collect(val))))
    end
end

labels2(df) = Dict(Pair.(propertynames(df), labels(df)))
