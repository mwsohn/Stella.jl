# function memory_saved(dt,n)

#     # total size and levels
#     tsize = sum([ismissing(x) ? 0 : length(x.value) for x in dt.data[n]])
#     s = Set(dt.data[n][:].values)

#     # for refs size
#     refsize = length(s) < typemax(UInt8) ? 1 : length(s) < typemax(UInt16) ? 2 : 4

#     # approximate size in categorical arrays
#     rsize = sum(length.(collect(s))) + refsize*dt.rows

#     return tsize*.8 > rsize ? true : false
# end

function get_numbytes(typelist,nvar)
    nb = Vector{UInt16}(undef,nvar)
    for i in 1:nvar
        if 0 < typelist[i] < 2045
            nb[i] = typelist[i]
        elseif typelist[i] == 32768
            nb[i] = 8 # 8 bytes
        elseif typelist[i] == 65526
            nb[i] = 8 # double
        elseif typelist[i] == 65527
            nb[i] = 4 # float
        elseif typelist[i] == 65528
            nb[i] = 4 # long
        elseif typelist[i] == 65529
            nb[i] = 2 # int
        elseif typelist[i] == 65530
            nb[i] = 1 # byte
        end
    end
    return nb
end

function strtonull(str::String)
	n = findfirst('\0',str)
	n == nothing && return str
	return str[1:n-1]
end

function alloc_array(vtype,vfmt,nobs::Int64)

    # create an Array for the relevant type
    if 0 <= vtype < 2045 || vtype == 32768 # string variable
        return Vector{Union{Missing,String}}(missing,nobs)
    elseif vtype == 65526
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(missing,nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return Vector{Union{Missing,DateTime}}(missing,nobs)
        else
            return Vector{Union{Missing,Float64}}(missing,nobs)
        end
    elseif vtype == 65527
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(missing,nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return Vector{Union{Missing,DateTime}}(missing,nobs)
        else
            return Vector{Union{Missing,Float32}}(missing,nobs)
        end
    elseif vtype == 65528
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(missing,nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return Vector{Union{Missing,DateTime}}(missing,nobs)
        else
            return Vector{Union{Missing,Int32}}(missing,nobs)
        end
    elseif vtype == 65529
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(missing,nobs)
        else
            return Vector{Union{Missing,Int16}}(missing,nobs)
        end
    elseif vtype == 65530
        return Vector{Union{Missing,Int8}}(missing,nobs)
    end

    error(vtype, " is not a valid variable type in Stata.")
end

"""
	read_stata(fn::String; chunks::Int=10, read_labels=false)

converts a stata datafile `fn` to Julia DataFrame. The original data file bigger than 100MB will be read in `chunks` (default = 10)
to save memory. If `read_labels` is set to `true`, it will not convert the data but will extract labels (both variable labels
and value labels) from the stata data file and convert them to Julia [Lables](https://github.com/mwsohn/Labels.jl).
"""
function read_stata(fn::String; chunks::Int=10, read_labels=false)

    fh = open(fn,"r")

    # dta file (<stata_dta><header><release>)
    header = String(read(fh,67))
    if header[2:10] != "stata_dta"
        error("Not a version 13 or 14 data file")
    end

    # data format version
    release = parse(Int16,header[29:31])
	if release == 117 # version 13
		len_varname = 33
		len_format = 49
		len_labelname = 33
		len_varlabel = 81
	elseif release == 118 # version 14
		len_varname = 129
		len_format = 57
		len_labelname = 129
		len_varlabel = 321
	else
		error("Can't convert data format version ",release,".")
	end

    # byte order: LSF or MSF
    byteorder = header[53:55]
    if byteorder == "MSF"
        error("Big-endian data are not supported yet.")
    end

    # number of variables
    skip(fh,3) # <K>
    nvar = read(fh,Int16)

    # number of observations
    skip(fh,7) #</K><N>
    nobs = Int(read(fh,Int32))

    # dataset label length
	if release == 117
		skip(fh,11)
		dslabel_len = read(fh,Int8)
    elseif release == 118
		skip(fh,15)
		dslabel_len = read(fh,Int16)
	end

    # read the label
    if dslabel_len > 0
        dslabel = String(read(fh,dslabel_len))
    end

    # time stamp
    skip(fh,19)
    timestamplen = read(fh,Int8)
    if timestamplen > 0
        timestamp = String(read(fh,timestamplen))
    end

    # map
    skip(fh,26) # </timestamp></header><map>
    statamap =Vector{Int64}(undef,14)
    read!(fh,statamap)

    # variable types
    skip(fh,22) # </map><variable_types>
    typelist = Vector{UInt16}(undef,nvar)
    read!(fh,typelist)

    # variable names
    skip(fh,27)
    varlist = Vector{Symbol}(undef,nvar)
    for i in 1:nvar
        varlist[i] = Symbol(strtonull(String(read(fh,len_varname))))
    end

    # sort list
    skip(fh,21) # </varnames><sortlist>
    for i in 1:nvar
        srt = read(fh,Int16)
    end

    # formats
    skip(fh,22) # </sortlist><formats> + 2 (2 bytes left over from the previous sequence)
    fmtlist = Vector{String}(undef,nvar)
    for i in 1:nvar
        fmtlist[i] = strtonull(String(read(fh,len_format)))
    end

    # value label names
    skip(fh,29) # </formats><value_label_names>
    valuelabels = Vector{String}(undef,nvar)
    numvlabels = 0
    for i in 1:nvar
        valuelabels[i] = strtonull(String(read(fh,len_labelname)))

        # count the number of value labels
        if length(valuelabels[i]) > 0
            numvlabels += 1
        end
    end

    # variable labels
    skip(fh,37) # </value_label_names><variable_labels>
    varlabels = Vector{String}(undef,nvar)
    for i in 1:nvar
        varlabels[i] = strtonull(String(read(fh,len_varlabel)))
    end

    # characteristics - we will not import them
    # read until we hit '</characteristics>'
    skip(fh,35) # </variable_labels><characteristics>
    while (true)
        readuntil(fh,'<')
        if String(copy(read(fh,5))) == "data>"
            break
        end
    end

    # nubmer of bytes for each variable
    numbytes = get_numbytes(typelist,nvar)

    # total length of each observation
    rlen = sum(numbytes)

    # number of bytes to skip in IOBuffer
    numskip = zeros(Int,nvar)
    numskip[1] = 0
    for i in 2:length(numbytes)
        numskip[i] = numbytes[i-1] + numskip[i-1]
    end

	# save the start position of the data section
	data_pos = position(fh)

	# skip the data section for now
	skip(fh,rlen*nobs)

    # if there is strLs, read them now
    skip(fh,7) # </data>
    tst = String(read(fh,7))
    if tst == "<strls>"
        # read strLs
        #
        # strL. Stata 13 introduced long strings up to 2 billon characters. strLs are
        # separated by "GSO".
        # (v,o): Position in the data.frame.
        # t:     129/130 defines whether or not the strL is stored with a binary 0.
        # len:   length of the strL.
        # strl:  long string.
        strls = OrderedDict()
        t = OrderedDict()
        while (String(read(fh,3)) == "GSO")
            v = read(fh,Int32)
            if release == 117
                o = read(fh,Int32)
            elseif release == 118
                o = read(fh,Int64)
            end
            t[(v,o)] = read(fh,UInt8)
            len = read(fh,UInt32)
            strls[(v,o)] = String(read(fh,len))
        end
    end

	tst = String(read(fh,5))
    if tst == "trls>"
        skip(fh,14)
    else
        error("Wrong position")
    end

	if read_labels == true

		# read value labels
    	# loop through all value labels
    	# define a Dict first

    	value_labels = Dict()

    	# numvlabels is defined above at the header section
    	for i in 1:numvlabels

			skipstr = String(read(fh,5))
			if skipstr != "<lbl>"
				break
			end
        	len = read(fh,Int32)
        	labname = Symbol(strtonull(String(read(fh,len_labelname))))

        	skip(fh,3) # padding
	        numvalues = read(fh,Int32) # number of entries
	        txtlen = read(fh,Int32) # length of value label text
	        value_labels[labname] = Dict()

	        #
	        offset = Vector{Int32}(undef,numvalues)
			read!(fh,offset) # offset
			values = Vector{Int32}(undef,numvalues)
	        read!(fh,values) # values
	        valtext = String(read(fh,txtlen)) # text table

	        for k in 1:numvalues
	            if k == numvalues
	                offset_end = txtlen
	            else
	                offset_end = offset[k+1]
	            end
	            if values[k] != ""
	                value_labels[labname][values[k]] = strtonull(valtext[offset[k]+1:offset_end])
	            end
	        end
	        skip(fh,6) # </lbl>
	    end

		variable_dict = Dict()
		lblname_dict = Dict()
	    for i in 1:nvar
	        variable_dict[varlist[i]] = varlabels[i]
	        lblname_dict[varlist[i]] = Symbol(valuelabels[i])
	    end

		return Label(variable_dict,value_labels,lblname_dict)
	end

	# read data now
	seek(fh,data_pos)

	# if the data size < 100MB, then
	# slurp the entire data section into memory
	# otherwise, we will read the data by smaller batches
	io = IOBuffer()
	if rlen*nobs < 100_000_000
		write(io,read(fh,rlen*nobs))
		seek(io,0)
		rdf = _read_dta(io,release,rlen,nobs,nvar,varlist,typelist,fmtlist,numskip,strls)
	else
		len = max(100000,ceil(Int,nobs/chunks))
		totlen = nobs
		rdf = DataFrame()
		for nread in 1:ceil(Int,nobs/len)
			if len < totlen
				totlen -= len
			else
				len = totlen
			end
			seek(io,0)
			write(io,read(fh,rlen*len))
			seek(io,0)

			rdf = vcat(rdf,_read_dta(io,release,rlen,len,nvar,varlist,typelist,fmtlist,numskip,strls))
		end
	end

	# close the file
	close(fh)

    # attach labels
    # data label
    if length(dslabel) > 0
        Stella.set_data_label!(rdf,dslabel)
    end

    # variable labels
    for i in 1:nvar
        Stella.set_col_label!(rdf,Symbol(varlist[i]),varlabels[i])
    end

    # value labels

	return rdf
end

function read_labels(fn::String)
	return read_stata(fn,read_labels=true)
end

function _read_dta(io, release, rlen, len, nvar,varlist,typelist,fmtlist,numskip,strls)

	df = DataFrame()

	dataitemf32::Float32 = 0.
	dataitemf64::Float64 = 0.
	dataitemi8::Int8 = 0
	dataitemi16::Int16 = 0
	dataitemi32::Int32 = 0
	v::Int32 = 0
	o::Int64 = 0
	z = zeros(UInt8,8)

	# interate over the number of variables
	for j in 1:nvar

		df[!,varlist[j]] = alloc_array(typelist[j],fmtlist[j],len)

		# len == number of observation in the batch
		for i in 1:len

			seek(io,numskip[j] + (i-1)*rlen)

			if 0 <= typelist[j] < 2045
				df[i,j] = strtonull(String(read(io,typelist[j])))
				# if empty string, return missing
				if df[i,j] == ""
					df[i,j] = missing
				end
			elseif typelist[j] == 32768 # long string
				if release == 117
					v = read(io,Int32)
					o = read(io,Int32)
				elseif release == 118
					z = read(io,8)
					v = reinterpret(Int16,z[1:2])[1]
					o = (reinterpret(Int64,z)[1] >> 16)
				end
				if (v,o) == (0,0)
					df[i,j] = missing
				else
					df[i,j] = strtonull(strls[(v,o)])
				end
			elseif typelist[j] == 65526
				dataitemf64 = read(io,Float64)
				if dataitemf64 > 8.9884656743e307
					df[i,j] = missing
				elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
					# convert it to Julia date
					df[i,j] = Date(1960,1,1) + Dates.Day(round(Int,dataitemf64))
				elseif fmtlist[j][1:3] == "%tc" || fmtlist[j][1:3] == "%tC"
					df[i,j] = DateTime(1960,1,1,0,0,0) + Dates.Millisecond(round(Int,dataitemf64))
				else
					df[i,j] = dataitemf64
				end
			elseif typelist[j] == 65527
				dataitemf32 = read(io,Float32)
				if dataitemf32 > 1.70141173319e38
					df[i,j] = missing
				elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
					# convert it to Julia date
					df[i,j] = Date(1960,1,1) + Dates.Day(round(Int,dataitemf32))
				elseif fmtlist[j][1:3] == "%tc" || fmtlist[j][1:3] == "%tC"
					df[i,j] = DateTime(1960,1,1,0,0,0) + Dates.Millisecond(round(Int,dataitemf32))
				else
					df[i,j] = dataitemf32
				end
			elseif typelist[j] == 65528
				dataitemi32 = read(io,Int32)
				if dataitemi32 > 2147483620
					df[i,j] = missing
				elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
					# convert it to Julia date
					df[i,j] = Date(1960,1,1) + Dates.Day(dataitemi32)
				elseif fmtlist[j][1:3] == "%tc" || fmtlist[j][1:3] == "%tC"
					df[i,j] = DateTime(1960,1,1,0,0,0) + Dates.Millisecond(dataitemi32)
				else
					df[i,j] = dataitemi32
				end
			elseif typelist[j] == 65529
				dataitemi16 = read(io,Int16)
				if dataitemi16 > 32740
					df[i,j] = missing
				elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
					# convert it to Julia date
					df[i,j] = Date(1960,1,1) + Dates.Day(dataitemi16)
				else
					df[i,j] = dataitemi16
				end
			elseif typelist[j] == 65530
				dataitemi8 = read(io,Int8)
				if dataitemi8 > 100
					df[i,j] = missing
				else
					df[i,j] = dataitemi8
				end
			end
		end
		# strls will be converted to categorical regardless of `categorize` option
		if typelist[j] == 32768
			categorical!(df,varlist[j])
		end

        # for vectors without missing values, convert the vector to an appropirate type
        if sum(ismissing.(df[!,varlist[j]])) == 0
            df[!,varlist[j]] = convert(Vector{nonmissingtype(eltype(df[!,varlist[j]]))},df[!,varlist[j]])
        end
	end

	return Stella.dfcompress(df)
end

# """
#     df, label = read_stata(fn::String)

# Converts a Stata datafile into a Julia DataFrame. It produces two memory objects (`df` and `label`).
# The first is the dataframe. The second is a [Lables](https://github.com/mwsohn/Labels.jl) object
# that provide functionality to attach variable and value labels.
# """
# function read_stata(fn::String)

#     dt = ReadStat.read_dta(fn)
#     df = readstat2dataframe(dt)

#     return df,get_labels(dt)
# end


#############################################################################
#
# DataFrame tools
#
#############################################################################
"""
    dfcompress(df::DataFrame)

Reduce `df`'s memory use by changing the eltype of each column in the `df` to
the type that can accommodate the larget and the smallest values within the same
integer or float class.
"""
function dfcompress(df::DataFrame; verbose = false)

    df2 = deepcopy(df)

    for v in propertynames(df)

        # if Array is empty after all missings are dropped
        # drop it from the df
        if all(ismissing.(df2[!,v]))
            if verbose
                select!(df2,Not(v))
                println(v, " was empty, now deleted")
            end                
            continue
        end

        # get the original eltype
        eltype_old = nonmissingtype(eltype(df2[!,v]))

        # if string, continue
        if eltype_old == String
            continue
        end

        # compress
        df2[!,v] = Stella.acompress(df2[!,v])

        if verbose && eltype_old != nonmissingtype(eltype(df2[!,v]))
            println(v, " was ", eltype_old, ", now ", nonmissingtype(eltype(df2[!,v])))
        end
    end

    return df2
end

"""
    acompress(da::AbstractVector)

Compresses a vector to a smallest numeric type that can hold without loss of information.
"""
function acompress(da::AbstractVector)

    # get the original eltype
    eltyp = eltype(da)
    eltype_old = nonmissingtype(eltyp)
    
    nomiss = true
    if eltyp == Union{Missing,eltype_old} && sum(ismissing.(da)) > 0
        nomiss = false;
    end

    # string variable - do not compress
    if eltype_old == String
        return da
    end

    if  eltype_old <: Integer

        varmin = minimum(collect(skipmissing(da)))
        varmax = maximum(collect(skipmissing(da)))
        if eltype_old != Int8 && varmin <= typemin(Int8) && varmax <= typemax(Int8)
            if nomiss
                return convert(Vector{Int8},da)
            else
                return convert(Vector{Union{Missing,Int8}},da)
            end
        elseif eltype_old != Int16 && varmin <= typemin(Int16) && varmax <= typemax(Int16)
            if nomiss
                return convert(Vector{Int16}, da)
            else
                return convert(Vector{Union{Missing,Int16}}, da)
            end
        elseif eltype_old != Int32 && varmin <= typemin(Int32) && varmax <= typemax(Int32)
            if nomiss
                return convert(Vector{Int32}, da)
            else
                return convert(Vector{Union{Missing,Int32}}, da)
            end
        end
    elseif eltype_old <: AbstractFloat
        # first test if the floats are integer numbers
        if all(collect(skipmissing(da)) .% 1 == 0) == true
            if nomiss
                return acompress(convert(Vector{Int64},da))
            else
                return acompress(convert(Vector{Union{Missing,Int64}},da))
            end
        end
    end
    return da
end


function atype(df::DataFrame,v::Symbol)
    # Array type = DA for DataArray, CA for Categorical Array, and UV for Union Vector
    if isdefined(Main,:CategoricalArrays) && typeof(df[!,v]) <: CategoricalArray
        return string("CA (", replace(string(eltype(df[!,v].refs)),"UInt" => ""), ")")
    # elseif isdefined(Main,:DataArrays) && typeof(df[!,v]) <: DataArray
    #      return "DA"
    elseif isdefined(Main,:PooledArrays) && typeof(df[!,v]) <: PooledArray
         return string("PA (", replace(string(eltype(df[!,v].refs)),"UInt" => ""), ")")
    # elseif isa(eltype(df[!,v]),Union)
    #     return "UV" # Union Vector
    end
    return ""
end

function etype(df::DataFrame,v::Symbol)
    # Eltype
    if typeof(df[!,v]) <: CategoricalArray
        eltyp = string(eltype(df[!,v].pool.levels))
        if in(eltyp,["String","AbstractString"])
            eltyp = string("Str",getmaxwidth(df[!,v].pool.levels))
        end
    else
        eltyp = string(nonmissingtype(eltype(df[!,v])))
        if in(eltyp,["String","AbstractString"])
            eltyp = string("Str",getmaxwidth(df[!,v]))
        elseif eltyp == "Dates.Date"
            eltyp = "Date"
        elseif eltyp == "Dates.DateTime"
            eltyp = "DateTime"
        end
    end

    return eltyp
end

function eltype2(df::DataFrame,v::Symbol)
    if typeof(df[!,v]) <: CategoricalArray
        return eltype(df[!,v].pool.index)
    end
    return nonmissingtype(eltype(df[!,v]))
end

"""
This is a section for attaching labels to a DataFrame using metadata and colmetadata.
The following functions will be named with "data_label", "col_label", 
"""

"""
    set_data_label(df::AbstractDataFrame,label::String)

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
    metadata(_df,"description")
    return nothing
end

"""
    delete_data_label(df::AbstractDataFrame)

Deletes the data label from the `df` DataFrame
"""
function delete_data_label!(_df::AbstractDataFrame)
    if "description" in metadatakeys(_df)
        deletemetadata!(_df,"description")
    end
    return nothing
end

"""
    set_fmtlib(df::AbstractDataFrame,path::String)

Saves the path to the format library file. A format file contains a dictionary whose keys are
format names in Symbols and values are a dictionary of value-description pairs. For example,
`race` variable has the following values: 1 for non-Hispanic Whites, 2 for non-Hispanic Blacks,
3 for Hispanics, and 4 for Other. Then one entry in this format libary may be:

fmtlib = Dict(
    :racelab = Dict(1 => "NH Whites", 2 => "NH Blacks", 3 => "Hispanics", 4 => "Other"),
    :agecat = Dict(1 => "< 40y", 2 => "40 - 64y", 3 => "65 - 74y", 4 => "75y or older")
)

Save this dictionary into a JLD2 file using

using JLD2
JLD2.save_object("fmtlib.jld2", fmtlib)

The location as a full or relative path for this file should be used as
the value for the `path` argument.
"""
function set_fmtlib!(_df::AbstractDataFrame,path::String)
    metadata!(_df,"fmtlib",path)
    return nothing
end

"""
    delete_fmtlib(df::AbstractDataFrame)

Deletes the fmtlib path from the `df` DataFrame.
"""
function delete_fmtlib!(_df::AbstractDataFrame)
    if "fmtlib" in metadatakeys(_df)
        deletemetadata!(_df,"fmtlib")
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
end

function delete_col_label!(_df::AbstractDataFrame,varname::Union{String,Symbol})
    if "label" in colmetadatakeys(_df,varname)
        deletecolmetadata!(_df,varname,"label")
    end
    return nothing
end



"""
    desc(df::DataFrame,varnames::Symbol...; labels::Union{Nothing,Label}=nothing, dfout::Bool = false, nmiss::Bool = false)

Displays variables in a dataframe much like `showcols`. It can display additional
attributes such as variable labels, value labels and display formats (not used in Julia)
if an optional `labels` is specified. It mimics Stata's `describe` command.
`labels` is automatically converted from a stata file by `read_stata` function. Or one can
be easily created as described in [Labels](https://github.com/mwsohn/Labels.jl).
"""
function desc(df::DataFrame,varnames::Symbol...; labels::Union{Nothing,Label}=nothing, dfout::Bool = false, nmiss::Bool = true)

    if length(varnames) == 0
        varnames = propertynames(df)
    end

    # number of variables
    numvar = length(varnames)

    # output dataframe
    dfv = DataFrame(Variable = varnames)
    dfv[!,:ArrayType] = Vector{String}(undef,size(dfv,1))
    dfv[!,:Eltype] = Vector{String}(undef,size(dfv,1))
    if nmiss
    	dfv[!,:Missing] = Vector{String}(undef,size(dfv,1))
    end

    if length(colmetadatakeys(df) ) > 0
        dfv[!,:Lblname] = Vector{String}(undef,size(dfv,1))
        dfv[!,:Description] = Vector{String}(undef,size(dfv,1))
    end

    for (i,v) in enumerate(collect(varnames))

        # variable name
        varstr = string(v)

	    # Array Type
        dfv[i,:ArrayType] = atype(df,v)

        # Eltype
        dfv[i,:Eltype] = etype(df,v)

        # percent missing
        if nmiss
            _nmiss = nmissing(df[!,v])
                dfv[i,:Missing] = string(round(100 * _nmiss/size(df,1),digits=1),"%")
        end

        if labels != nothing
            dfv[i,:Lblname] = lblname(labels,v) == nothing ? "" : string(lblname(labels,v))
        end

        varlabel = col_label(df)
        if length(varlabel) > 0
            dfv[i,:Description] = varlabel[v]
        end
    end
					
    header = ["Variable","Atype","Eltype"]
    alignment = [:l,:l,:l]
    if nmiss
	header = vcat(header,"% Miss")
	alignment = vcat(alignment,:r)
    end
    if length(varlabel) > 0
    	header = vcat(header,["Lbl Name","Description"])
	    alignment = vcat(alignment,[:l,:l])
    end

    if dfout 
    	return dfv
    else
        if "description" in metadatakeys(df)
            println(data_label(df))
        end
        pretty_table(dfv,
            alignment=alignment,
            header=header,
            crop=:none,
            show_row_number = true)
    end
end

function nmissing(s::AbstractArray)
    return sum(ismissing.(s))
end

function getmaxwidth(s::AbstractArray)
    if isa(s, CategoricalArray) && nonmissingtype(eltype(s)) <: CategoricalString
	return maximum(length.(s.pool.levels))
    end
	
    if nmissing(s) == size(s,1)
	return 0
    end
	
    return  maximum(length.(collect(skipmissing(s))))
end

"""
    ds(df::DataFrame, typ::Type, args...)
    ds(df::DataFrame, re::Regex)

Returns an array of Symbols of the columns whose eltype matches the `typ` or the regular expression `re`.
For `String` variables, up to two additional arguments can be specified to indicate the minimum and maximum widths
to be returned. In Julia, string variables does not have an associated widths as in Stata.
So, in this function, the width is measured by the length of the longest string. The function follows `Number` type
structure. For example, to list all numeric variables, use `Number`; for all integer variables, use `Integer`; for
all floating point numbers, use `AbstractFloat`. This function is useful to subset a DataFrame using
eltypes or column names. This function emulates the Stata `ds` command.

# Example
```
julia> showcols(aapl)
6081×7 DataFrames.DataFrame
│ Col # │ Name      │ Eltype  │ Missing │
├───────┼───────────┼─────────┼─────────┤
│ 1     │ Date      │ String  │ 0       │
│ 2     │ Open      │ Float64 │ 0       │
│ 3     │ High      │ Float64 │ 0       │
│ 4     │ Low       │ Float64 │ 0       │
│ 5     │ Close     │ Float64 │ 0       │
│ 6     │ Volume    │ Int64   │ 0       │
│ 7     │ Adj_Close │ Float64 │ 0       │

julia> ds(aapl,Int64)
1-element Array{Symbol,1}:
 :Volume

julia> ds(aapl,Float64)
5-element Array{Symbol,1}:
 :Open
 :High
 :Low
 :Close
 :Adj_Close

julia> ds(aapl,Number)
6-element Array{Symbol,1}:
 :Open
 :High
 :Low
 :Close
 :Volume
 :Adj_Close

julia> ds(aapl,String)
1-element Array{Symbol,1}:
 :Date

julia> ds(aapl,r".*lose")
2-element Array{Symbol,1}:
 :Close
 :Adj_Close

```
"""
function ds(df::DataFrame, typ::Type, args...)

    dslist = Vector{Symbol}()

    for v in propertynames(df)
        if typ == String && nonmissingtype(eltype(df[!,v])) == String
            if length(args) == 0
                push!(dslist,v)
                continue
            end

            maxlen = getmaxwidth(df[!,v])
            if length(args) == 1
                if args[1] <= maxlen
                    push!(dslist,v)
                end
            elseif length(args) == 2
                if args[1] <= maxlen <= args[2]
                    push!(dslist,v)
                end
            end
        elseif nonmissingtype(eltype(df[!,v])) == typ
            push!(dslist,v)
        elseif (typ == Integer && nonmissingtype(eltype(df[!,v])) <: Integer) || (typ == AbstractFloat && nonmissingtype(eltype(df[!,v])) <: AbstractFloat) || (typ == Number && nonmissingtype(eltype(df[!,v])) <: Number)
            push!(dslist,v)
        end
    end

    return dslist

end
function ds(df::DataFrame,re::Regex)
    dslist = Array{Symbol,1}()

    for v in propertynames(df)
        if occursin(re,string(v))
            push!(dslist,v)
        end
    end

    return dslist
end


"""
    pickone(df::DataFrame,groupvars::Vector{Symbol})

Create a vector that identifies one record in a `groupvars` in the `df` DataFrame.
This function creates a variable similar to the Stata command
`egen byte pickone = tag(groupvar)`. `groupvars` can be an array of Symbols or a single
Symbol.
"""
function pickone(df::DataFrame,groupvars::Vector{Symbol})
    df[:_____obs_____] = collect(1:size(df,1))
    done = zeros(Int8,size(df,1))
    for subdf in groupby(df, groupvars)
        done[subdf[1,:_____obs_____]]=1
    end
    select!(df,Not(:_____obs_____))
    return done
end
pickone(df::DataFrame,groupvar::Symbol) = pickone(df, [groupvar])


"""
    duplicates(df::DataFrame, args::Symbol...; cmd::Symbol = :report)

Reports, tags, or deletes duplicates in a DataFrame for one or more columns. Use `cmd` to
request one of the actions:

- :report - for getting a frequency table for the number of duplicates
- :tag - for identifying rows with duplicate values
- :drop - for dropping all rows with duplicate values except for the first row.
"""
function duplicates(df::DataFrame, args::Union{Symbol,String}... ; cmd::Symbol = :report)

    if in(cmd,[:report,:drop,:tag]) == false
        error("`cmd = `", cmd, "` is not supported.")
    end

    # if args are not specified, use all variables
    nargs = length(args)
    if nargs == 0
        args = tuple(names(df)...)
    end

    gdf = groupby(df,collect(args))
    dfx = combine(gdf,nrow => :__dups__)

    if cmd == :report
        na = freqtable(dfx, :__dups__)
        setdimnames!(na,"copies",1)
        return na
    end

    # substract 1 from :x in dfx (we are reporting 0 for unique observations)
    dfx[:__dups__] .-= 1
    df = leftjoin(df, dfx, on = collect(args))

    if cmd == :tag
        return df[:__dups__]
    end

    # if cmd == :drop
    ba = [x == 1 ? true : false for x in pickone(df,collect(args))]
    return select(df[ba,:],Not(:__dups__))
end

"""
    sample(df::AbstractDataFrame,num::Real)

Creates an Int8 array that identifies a randomly selected sample (1 = sample; 0 otherwise) from the input dataframe. `num`
specifies the sample size. If `num` is an integer, it indicates the number of rows to be selected.
If it is a float, it indicates the percentage of the input dataframe to be randomly selected.
Use `Random.seed!()` to set a seed before using `sample()`.

##Example
To select 100 rows, use

```
julia> df[!,:sample] = sample(df,100)
```

To select 20% of the original dataframe, use

```
julia> df[!,:sample2] = sample(df,.2)
```

"""
function sample(df::AbstractDataFrame,n::Real)
    len = size(df,1)
    asamp = zeros(Int8,len)
    if isa(n,AbstractFloat)
        asamp[randsubseq(collect(1:len),n)] .= 1
    else
        asamp[randperm(len)[1:n]] .= 1
    end
    return asamp
end

"""
    dfsample(df::AbstractDataFrame,num::Real)

Creates a dataframe with a randomly selected sample from the input dataframe. `num`
specifies the sample size. If `num` is an integer, it indicates the number of rows to be selected.
If it is a float, it indicates the percentage of the input dataframe to be randomly selected.
Use `Random.seed!()` to set a seed before running `dfsample()`.

##Example
To select 100 rows, use

```
julia> df2 = dfsample(df,100)
```

To select 20% of the original dataframe, use

```
julia> df2 = dfsample(df,.2)
```

"""
function dfsample(df::AbstractDataFrame,num::Real)
    return df[sample(df,num) .== 1,:]
end

"""
    dfmerge(df1::DataFrame,df2::DataFrame,linkers::Union{Symbol,Vector}; kind::Symbol = :outer)

Produces a merged dataframe with a `:_merge` variable indicating how the merge was done.
This function is a wrapper for `join` in the DataFrames package and mimics `merge ... using ...`
command in Stata. Currently, it supports a merge of two sources only.

The default merge is a `outer` join that keeps records from both sources in the merged
dataframe. `:_merge` values indicates:

- `1`: in the first source alone
- `2`: in the second source alone
- `3`: in both sources. A successful merge was achieved.

It is the responsibility of the user to keep or discard records with `:_merge` values 1 or 2.

The third argument can be either a single symbol or a vector of symbols to be used as the
linkage variables.

Optionally, `kind` argument can be specified. Available options are
`:left`, `:right`, `:inner`,`:outer`,`:semi`, and `:anti`. Consult the documentation for
the `join` function in the DataFrames package.
"""
function dfmerge(df1::DataFrame,df2::DataFrame,linkers::Union{Symbol,Vector})
    if in(:_merge,propertynames(df1))
        error("`:_merge' exists in the first dataframe")
    end
    if in(:_merge,propertynames(df2))
        error("`:_merge' exists in the second dataframe")
    end

    df1[!,:___mergeleft___] = ones(Int8,size(df1,1))
    df2[!,:___mergeright___] = ones(Int8,size(df2,1))

    df_merged = outerjoin(df1,df2,on = linkers)
    df_merged[!,:_merge] = zeros(Int8,size(df_merged,1))
    for i = 1:size(df_merged,1)
        if ismissing(df_merged[i,:___mergeright___])
            df_merged[i,:_merge] = 1
        elseif ismissing(df_merged[i,:___mergeleft___])
            df_merged[i,:_merge] = 2
        elseif df_merged[i,:___mergeleft___] == 1 && df_merged[i,:___mergeright___] == 1
            df_merged[i,:_merge] = 3
        end
    end
    select!(df_merged,Not([:___mergeleft___,:___mergeright___]))
    return df_merged
end

function get_class(val::Real,thresholds::Vector,lower::Bool)
    if lower == false
        for i = 1:length(thresholds)
            if i == 1 && val < thresholds[1]
                return 1
            elseif i > 1 && thresholds[i-1] <= val < thresholds[i]
                return i
            end
        end
    else
        for i = 1:length(thresholds)
            if i == 1 && val <= thresholds[1]
                return 1
            elseif i > 1 && thresholds[i-1] < val <= thresholds[i]
                return i
            end
        end
    end
    return length(thresholds)+1
end

"""
    classify(da::AbstractVector,thresholds::Vector; lower::Bool = false)

Produces a categorical variable based on a data array and a vector of thresholds.
Use `lower = true` to classify the threshold value to the lower class.
"""
function classify(da::AbstractVector,thresholds::Vector; lower::Bool = false)
    da2 = Vector{Union{Missing,Int8}}(undef,length(da))
    if length(thresholds) > 100
        error("Cannot use more than 100 threshold values.")
    end

    for i = 1:size(da,1)
        if ismissing(da[i])
            da2[i] = missing
        else
            da2[i] = get_class(da[i],thresholds,lower)
        end
    end
    return da2
end

"""
    addterms(fm::Formula,v::Vector{Symbol})

Adds covariates to an existing formula. `fm` is an object of `Formula` type
created by `@formula` macro. `v` is a vector of Symbols. This function is useful to
build hierarchically nested models.

### Example

```
julia> using DataFrames

julia> fm = @formula(income ~ age + male)
Formula: income ~ age + male

julia> fm2 = addterms(fm,[:race, :educ])
Formula: income ~ age + male + race + educ
```

"""
function addterms(fmm::Formula,v::Vector{Symbol})

    # create a new Formula object
    fm = deepcopy(fmm)

    # if fm.rhs is a null modle with rhs == 1 or
    # a unadjusted model with only one rhs variable,
    # convert it to an Expr first
    if fm.rhs == 1
        fm.rhs = Expr(:call,:+,Tuple(v)...)
    elseif typeof(fm.rhs) == Symbol
        fm.rhs = Expr(:call,:+,fm.rhs,Tuple(v)...)
    else
        push!(fm.rhs.args,Tuple(v)...)
    end

    return fm
end

"""
    Stella.formula(o::Symbol,v::Vector{Symbol})

Creates a Formula from a vector of covariates. This function is not exported for potential conflict with other packages.

- o - A symbol for the dependent variable
- v - A vector of symbols for covariates

## Usage

julia> fm = Stella.formula(:readmit, [:age,:male,:race,:chd,:exercise])
Formula: readmit ~ age + male + race + chd + exercise
"""
function formula(o::Symbol,v::Vector{Symbol})
    if length(v) == 1
        rhs = v[1]
    else
        rhs = Expr(:call,:+,Tuple(v)...)
    end
    e = :($o ~ $rhs)
    return Formula(e,e,o,rhs)
end

"""
    renvars!(df::DataFrame; vars::Array{Symbol,1}, case = "lower")

Rename column names in `vars` to either upper or lower cases. The default is to convert
all columns to lower case names.
"""
function renvars!(df::DataFrame; vars=[], case="lower")
    numvar = length(vars)
    symnames = propertynames(df)

    if numvar == 0
        varnames = propertynames(df)
    else
        varnames = propertynames(df[vars])
    end

    for nm in varnames

        if case == "lower" || case == "LOWER"
            newname = lowercase(string(nm))
        elseif case == "upper" || case == "UPPER"
            newname = uppercase(string(nm))
        else
            error(option," ", case," is not allowed.")
        end
        if nm != Symbol(newname)
            rename!(df,nm => Symbol(newname))
        end
    end
end

"""
	ci(df::AbstractDataFrame,colname::Union{Symbol,String})
	ci(df::AbstractDataFrame, colname1::Union{Symbol,String}, colname2::Union{Symbol,String})
	ci(df::AbstractDataFrame,colre::Regex)

Produces column index using DataFrames.columnindex() function. It can produce a single index,
or a range of indices, or multiple indices matching a regular expression.
"""
function ci(df::AbstractDataFrame,colname::Union{Symbol,String})
    return columnindex(df,colname)
end
function ci(df::AbstractDataFrame, colname1::Union{Symbol,String}, colname2::Union{Symbol,String})
    ci1 = columnindex(df,colname1)
    ci2 = columnindex(df,colname2)
    return ci1 < ci2 ? range(ci1,ci2) : range(ci2,ci1)
end
function ci(df::AbstractDataFrame,colre::Regex)

    vec = Int[]
    for n in names(df)
        m = match(colre, n)
        if m !== nothing
            push!(vec,columnindex(df,n))
        end
    end

    return vec
end


"""
    destring(da::AbstractArray;force=true)
    destring(df::DataFrmae,strvar::Symbol;force=true)
    destring!(df::DataFrame,strvars; newvars::Vector{Symbol} = [], force=true, replace=false)

Convert a string vector to a numeric vector. Use `force = true` to coerce conversion of alphanumeric strings to
missing values. If `force = false`, any vector containing non-numeric characters are not converted.
Use `replace = true` in `destring!` to replace the original string vector with a new converted numeric vector.
If `replace` option is specified, `newvars` array is ignored.
"""
function destring(da::AbstractArray; force=true)
    if length(da) == 0
        error("Input array is empty!")
    end
    if nonmissingtype(eltype(da)) <: Number
        error("Input array is a numeric vector. No conversion needed.")
    end

    # check if the values include any alphabetic characters or decimals
    isfloat = false
    alpha = false
    da_safe = collect(skipmissing(da))
    for i in length(da_safe)
        if sum([isdigit(x) == false for x in da_safe[i]]) > 0
            alpha = true
        end
        if occursin(r"[,0-9]*\.[0-9]+",da_safe[i])
            isfloat = true
        end
    end

    if alpha && force == false
        error("Input array contains alphabetic letters. Use 'force=true' option to coerce conversion.")
    end

    T = isfloat ? Float64 : Int64
    da2 = [ismissing(x) ? missing : parse(T,x) for x in da]

    return acompress(da2)
end
destring(df::AbstractDataFrame,strvar::Symbol; force=true) = destring(df[:,strvar],force=force)
function destring!(df::DataFrame,strvars; newvars::Vector{Symbol} = [], force=true, replace=false)

    if replace
        for v in strvars
            df[v] = destring(df[v],force=force)
        end
    else
        # check if there are same number of symbols in strvars
        if length(strvars) != length(newvars)
            error("The number of symbols in ", strvars, " and ", newvars, " are not the same.")
        end
        for i in 1:length(strvars)
            df[newvars[i]] = destring(df[strvars[i]], force=force)
        end
    end
end

# """
#     rowsum(df::DataFrame)
#
# Creates a DataArray that contains the row total of all values on the same row of `df`.
# If one of the columns contain an NA value on a row, an NA value will be returned for that
# row. This function emulates Stata's `egen rowsum = rowtotal(var1 - var3)`.
#
# ```
# julia>df[:rowsum] = df[[:var1,:var2,:var3]]
# ```
#
# If the position numbers for `:var1` (e.g., 4), `:var2` (5), `:var3` (6) are known and consecutive,
# you can specify them as follows:
#
# ```
# julia>df[:rowsum] = df[collect(4:6)]
# ```
# """
# function rowsum(df::DataFrame)
#
#     isfloat = false
#     for i in 1:size(df,2)
#         if eltype(df[i]) <: AbstractFloat
#             isfloat = true
#         end
#     end
#
#     if isfloat
#         da = zeros(Union{Missing,Float64},size(df,1))
#     else
#         da = zeros(Union{Missing,Int64},size(df,1))
#     end
#
#     ba = completecases(df)
#     for i = 1:size(df,1)
#         if ba[i] == false
#             da[i] = missing
#         else
#             for j = 1:size(df,2)
#                 da[i] += df[i,j]
#             end
#         end
#     end
#     return da
# end

"""
    rowstat(df::DataFrame,func::Function)

Creates a vector that contains the row statistic of all values on the same row of `df`
produced by the `func` function. If one of the columns contain a missing value on a row, the returned value will be mssing for that
row. This function emulates Stata's `egen` row functions such as `rowtotal`, `rowmean`, etc.

```
julia>df[:rowmean] = rowstat(df[[:var1,:var2,:var3]],mean)
```
"""
function rowstat(df::DataFrame,func::Function)

    # df size
    nrow, ncol = size(df)

    # array to be returned
    da = zeros(Union{Missing,Float64},size(df,1))

    # temporary array of row values
    ta = Vector{Union{Missing,Float64}}(undef,size(df,2))

    # iterate over all rows
    for i = 1:nrow

        # counter for non-missing values in each row
        k=0
        for j = 1:ncol
            if ismissing(df[i,j]) == false
                k += 1
                ta[k] = df[i,j]
            end
        end

        # if any of the columns has a missing value, return missing
        if k < ncol
            da[i] = missing
        else
            tmpfloat = func(ta[1:k])
            da[i] = isnan(tmpfloat) ? missing : tmpfloat
        end
    end
    return da
end

"""
    xtile(da::AbstractArray;nq::Int = 4, lower::Bool = false, cutoffs::Union{Nothing,AbstractVector} = nothing)
    xtile(df::DataFrame,varname::Symbol;nq::Int = 4, lower::Bool = false, cutoffs::Union{Nothing,AbstractVector} = nothing)

Create a DataArray{Int8,1} that identifies `nq` categories based on values in `da`.
The default `nq` is 4. `cutoffs` vector can be provided to make custom categories.
`cutoffs` vector is expected to contain `nq - 1` elements. The minimum and maximum values
will be computed. `lower = true` will make the threshold values to go in the lower category (default: false).

```
julia> df[:agecat] = xtile(df[:age], nq = 3)
```
"""
function xtile(da::AbstractArray ; nq::Int = 4, lower::Bool = false, cutoffs::Union{Nothing,AbstractVector} = nothing)

    if nq == 1
        error("`nq` must be greater than 1")
    end

    da2 = collect(skipmissing(da))
    if cutoffs == nothing
	    cutoffs = nquantile(da2,nq)
    elseif length(cutoffs) == nq - 1
        cutoffs = vcat(minimum(da2), cutoffs, maximum(da2))
    else
        error("`cutoffs` vector length is not consistent with `nq`. It must be 1 greater or 1 less than `nq`.")
    end

    if lower
        return [ismissing(x) ? missing : qval_low(x,cutoffs) for x in da2]
    else
        return [ismissing(x) ? missing : qval_high(x,cutoffs) for x in da2]
    end
end
xtile(df::DataFrame,arg::Symbol; nq::Int = 4, cutoffs::Union{Nothing,AbstractVector} = nothing) = xtile(df[!,arg], nq = nq, cutoffs = cutoffs)

function qval_high(val::Real,cut::Vector)
    cl = length(cut)

    for i in 2:cl-1
        if i == 2 && val < cut[i]
            return 1
        elseif i == cl-1 && cut[i] <= val
            return i
        elseif cut[i] <= val < cut[i+1]
            return i
        end
    end
    warn("Error - check qval function")
end

function qval_low(val::Real,cut::Vector)
    cl = length(cut)

    for i in 2:cl-1
        if i == 2 && val <= cut[i]
            return 1
        elseif i == cl-1 && cut[i] <= val
            return i
        elseif cut[i] < val <= cut[i+1]
            return i
        end
    end
    warn("Error - check qval function")
end

# """
#     rapply(df::DataFrame,r::OrderedDict)

# Creates a recoded vector of values according to the rule set specified as an ordered dictionary.
# Keys in the rule set must be the recoded value and the values must be the rules.
# The rule must be an expression or a string that return a boolean value (`true` or `false`).

# ## Example

# julia> ruleset = OrderedDict(
#             1 => :( df[:race] .== 1 && df[:hispanic] .== 0),
#             2 => :( df[:race] .== 2 && df[:hispanic] .== 0),
#             3 => :( df[:hispanic] .== 1),
#             4 => :( df[:hispanic] .== 0 && in.(df[:race],[1,2]) == false)
#        )
# OrderedDict{Int64,Expr} with 4 entries:
#   1 => :(df[:race] .== 1 && df[:hispanic] .== 0)
#   2 => :(df[:race] .== 2 && df[:hispanic] .== 0)
#   3 => :(df[:hispanic] .== 1)
#   4 => :(df[:hispanic] .== 0 && in.(df[:race], [1, 2]) == false)

# julia> rapply(adf,ruleset)


# """
# function rapply(df::DataFrame,r::OrderedDict)

#     # values
#     vals = sort(collect(keys(r)))

#     # types
#     vtyp = eltype(vals)

#     # empty Vector
#     vec = Vector{Union{vtyp,Missing}}(undef,size(df,1))

#     # go through the rules and assign values
#     for v in vals

#         if typeof(r[v]) == String
#             ba = eval(parse(r[v]))
#         elseif typeof(r[v]) == Expr
#             ba = eval(r[v])
#         else
#             error(typeof(r[v]), " is not an allowed type for rule ", string(r[v]))
#         end

#         for i in findall(x->x==true,ba)
#             vec[i] = v
#         end
#     end

#     return vec
# end



"""
	categorical!(df::AbstractDataFrame,vv::Union{Symbol,Vector{Symbol}})

converts vectors into CategoricalArrays
"""
function categorical!(df::AbstractDataFrame,vv::Union{Symbol,Vector{Symbol}})
    for v in vcat(vv)
        if isa(df[:,v],CategoricalArray) == false
            df[!,v] = categorical(df[!,v])
        end
    end
end


"""
	uncategorical!(df::AbstractDataFrame,vv::Union{Symbol,Vector{Symbol}})

uncategorizes a vector of CategoricalArrays by replacing the refs with their original values.
"""
function uncategorical!(df::AbstractDataFrame,vv::Union{Symbol,Vector{Symbol}})
    for v in vcat(vv)
        if isa(df[:,v],CategoricalArray)
            lev = df[:,v].pool.levels
            df[!,v] = [x == 0 ? missing : lev[x] for x in df[:,v].refs]
        end
    end
end

"""
	uncategorize(v::AbstractArray)

produces an array by uncategorizing a CategoricalArray by replacing the refs with its original values.
"""
function uncategorize(v::CategoricalArray)
    return [x == 0 ? missing : v.pool.levels[x] for x in v.refs]
end

# convert the Arrow data types from the feather file
# function convert_feather(df::DataFrame)
#     for v in propertynames(df)
#         if typeof(df[!,v]) <: DictEncoding
#             df[!,v] = CategoricalArray(df[!,v])
# 	elseif eltype(df[!,v]) == Feather.Arrow.Datestamp
# 	    df[!,v] = convert.(Date,df[!,v])		
#         else
#             df[!,v] = Array(df[!,v])
#         end
#     end

#     return df
# end


