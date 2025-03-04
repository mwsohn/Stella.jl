@reexport import TableMetadataTools

function labels!(df::AbstractDataFrame,dd::Dict)
    for key in keys(dd)
        label!(df,key,dd[key])
    end
end

function values!(df::AbstractDataFrame,v::Union{Symbol,String},val::Union{Dict,NamedTuple}; ordered = true)

    # if isa(df[:,v],CategoricalArray) 
    #     # uncategorical!(df,v)
    #     levels!(df[!, v], last.(sort(collect(val))))
    # elseif nonmissingtype(eltype(df[:,v])) <: Integer
    #     categorical!(df[!,v], ordered = ordered)
    #     recode!(df[:,v],val...)
    #     levels!(df[!,v], last.(sort(collect(val))))
    # end
    if isa(df[!, v], CategoricalArray)
        levels!(df[!, v], last.(sort(collect(val))))
    elseif nonmissingtype(eltype(df[:, v])) <: Integer
        df[!, v] = recode(df[:, v], val...)
        df[!, v] = categorical(df[!, v], ordered=ordered)
    end
end

labels2(df) = Dict(Pair.(propertynames(df), labels(df)))
