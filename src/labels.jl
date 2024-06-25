@reexport import TableMetadataTools

function labels!(df::DataFrame,dd::Dict)
    for key in keys(dd)
        label!(df,key,dd[key])
    end
end

function values!(df::DataFrame,v::Union{Symbol,String},valdict::Dict)
    if nonmissingtype(eltype(df[:,v])) <: Integer
        vv = categorical(df[!,v])
        df[!,v] = recode(vv,valdict...)
        levels!(df[!,v], last.(sort(collect(valdict))))
    end
end
