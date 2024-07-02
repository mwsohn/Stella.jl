@reexport import TableMetadataTools

function labels!(df::DataFrame,dd::Dict)
    for key in keys(dd)
        label!(df,key,dd[key])
    end
end

function values!(df::DataFrame,v::Union{Symbol,String},val::Union{Dict,NamedTuple}; ordered = true)
    if nonmissingtype(eltype(df[:,v])) <: Integer
        vv = categorical(df[!,v], ordered = ordered)
        df[!,v] = recode(vv,val...)
        levels!(df[!,v], last.(sort(collect(val))))
    end
end
