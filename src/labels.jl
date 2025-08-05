@reexport import TableMetadataTools

function labels!(df::AbstractDataFrame,dd::Dict)
    for key in keys(dd)
        label!(df,key,dd[key])
    end
end

function values!(df::AbstractDataFrame,v::Union{Symbol,String},val::Union{Dict,NamedTuple}; ordered = false)
    isa(df[!, v], CategoricalArray) && throw(ArgumentError("Already a Categorical Array. Use `levels!` to recode the level values"))
        
    if nonmissingtype(eltype(df[:, v])) <: Integer
        df[!, v] = recode(df[:, v], val...)
        df[!, v] = categorical(df[!, v], ordered=ordered)
    end
end

labels2(df) = Dict(Pair.(propertynames(df), labels(df)))
