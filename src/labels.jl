export set_data_label!, data_label, delete_data_label!,
    set_col_label!, col_label, delete_col_label!,
    set_lblname!, lblname, delete_lblname!,
    set_value_dict!, value_dict, delete_value_dict!,
    vallab

"""
    set_data_label!(df::AbstractDataFrame,label::String)

Saves a `label` to the `df` DataFrame. It does not return anything nor generates any
message.
"""
function set_data_label!(_df::AbstractDataFrame,label::String)
    metadata!(_df,"description",label,style=:default)
    return nothing
end

"""
    data_label(df::AbstractDataFrame)

Extracts the data label saved in the `df` DataFrame and returns it.
"""
function data_label(_df::AbstractDataFrame)
    if "description" in metadatakeys(_df)
        return metadata(_df,"description")
    end
    return nothing
end

"""
    delete_data_label!(df::AbstractDataFrame)

Deletes the data label from the `df` DataFrame
"""
function delete_data_label!(_df::AbstractDataFrame)
    if "description" in metadatakeys(_df)
        deletemetadata!(_df,"description")
    end
    return nothing
end

function set_value_dict!(_df::AbstractDataFrame,vlib::Dict)
    metadata!(_df,"value_dictionary",vlib)
    return nothing
end

function value_dict(_df::AbstractDataFrame)
    if "value_dictionary" in metadatakeys(_df)
        return metadata(_df,"value_dictionary")
    end
end

function delete_value_dict!(_df::AbstractDataFrame)
    if "value_dictionary" in metadatakeys(_df)
        deletemetadata(_df, "value_dictionary")
    end
    return nothing
end

function set_col_label!(_df::DataFrame,labels::Dict)
    colnames = names(_df)
    for (key,val) in labels
        if string(key) in colnames
            colmetadata!(_df,key,"label",val, style=:note)
        end
    end
    return nothing
end
function set_col_label!(_df::AbstractDataFrame,varname::Union{String,Symbol},label::String)
    colmetadata!(_df,varname,"label",label,style=:note)
    return nothing
end


function col_label(_df::AbstractDataFrame)
    label = Dict()
    for v in propertynames(_df)
        if "label" in colmetadatakeys(_df,v)
            label[v] = colmetadata(_df,v,"label")
        end
    end
    return label
end
function col_label(_df::AbstractDataFrame,varname::Union{String,Symbol})
    if "label" in colmetadatakeys(_df,varname)
        return colmetadata(_df,varname,"label")
    end
    return nothing
end

function delete_col_label!(_df::AbstractDataFrame,varname::Union{String,Symbol})
    if "label" in colmetadatakeys(_df,varname)
        deletecolmetadata!(_df,varname,"label")
    end
    return nothing
end

function set_lblname!(_df::AbstractDataFrame, varname::Union{Symbol,String}, vlabname::Symbol)
    colmetadata!(_df,varname,"lblname",vlabname)
    return nothing
end
function set_lblname!(_df::AbstractDataFrame,lblnames::Dict)
    for v in keys(lblnames)
        colmetadata!(_df, v, "lblname", lblnames[v])
    end
    return nothing
end

function lblname(_df::AbstractDataFrame, varname::Union{Symbol,String})
    if haskey(colmetadata(_df,varname),"lblname")
        return colmetadata(_df,varname,"lblname")
    end
end
function lblname(_df::AbstractDataFrame, varnames::AbstractVector)
    ldict = Dict()
    for v in varnames
        if haskey(colmetadata(_df,v),"lblname")
            ldict[v] = colmetadata(_df,v,"lblname")
        end
    end
    if length(ldict) > 0
        return ldict
    end
    return nothing
end
function lblname(_df::AbstractDataFrame)
    valdict = Dict()
    for v in propertynames(_df)
        if haskey(colmetadata(_df, v), "lblname")
            lblname = colmetadata(_df, v, "lblname")
            if lblname != ""
                valdict[v] = lblname
            end
        end
    end
    return valdict
end

function delete_lblname!(_df::AbstractDataFrame, varname::Union{Symbol,String})
    if haskey(colmetadata(_df, varname), "lblname")
        return deletecolmetadata(_df, varname, "lblname")
    end
    return nothing
end

function vallab(_df::AbstractDataFrame, varname::Union{Symbol,String})
    valdict = value_dict(_df)
    lname = lblname(_df,varname)
    if length(valdict) > 0 && haskey(valdict,lname)
        return Dict(varname => valdict[lname])
    end
    return nothing
end
function vallab(_df::AbstractDataFrame, varnames::AbstractVector)
    valdict = value_dict(_df) 
    lname = lblname(_df,varnames)
    if length(valdict) == 0 || length(lname) == 0
        return nothing
    end
    vdict = Dict(x => valdict[lname[x]] for x in keys(lname))
    if length(vdict) > 0
        return vdict
    end
    return nothing
end



