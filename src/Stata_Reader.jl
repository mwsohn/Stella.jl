using DataFrames, DataStructures, CategoricalArrays, Missings

function get_numbytes(typelist,nvar)
    nb = Array{UInt16,1}(nvar)
    for i in 1:nvar
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
end

function strtonull(str)

    if UInt32(str[1]) == 0x00
        return ""
    end

    for i in 1:endof(str)
        if isvalid(str,i) && UInt32(str[i]) == 0x00
            return str[1:i-1]
        end
    end

    return str
end

function alloc_dataarray(vtype,vfmt,nobs::Int64)

    # create a DataArray for the relevant type
    if 0 <= vtype < 2045 || vtype == 32768 # string variable
        return DataArray(String,nobs) #Vector{Union{Missing,String}}(nobs)
    elseif vtype == 65526
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return DataArray(Date,nobs) # Vector{Union{Missing,Date}}(nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return DataArray(DateTime,nobs) #Vector{Union{Missing,DateTime}}(nobs)
        else
            return DataArray(Float64,nobs) #Vector{Union{Missing,Float64}}(nobs)
        end
    elseif vtype == 65527
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return DataArray(Date,nobs) #Vector{Union{Missing,Date}}(nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return DataArray(DateTime,nobs) #Vector{Union{Missing,DateTime}}(nobs)
        else
            return DataArray(Float32,nobs) #Vector{Union{Missing,Float32}}(nobs)
        end
    elseif vtype == 65528
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return DataArray(Date,nobs) #Vector{Union{Missing,Date}}(nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return DataArray(DateTime,nobs) #Vector{Union{Missing,DateTime}}(nobs)
        else
            return DataArray(Int32,nobs) #Vector{Union{Missing,Int32}}(nobs)
        end
    elseif vtype == 65529
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return DataArray(Date,nobs) #Vector{Union{Missing,Date}}(nobs)
        else
            return DataArray(Int16,nobs) #Vector{Union{Missing,Int16}}(nobs)
        end
    elseif vtype == 65530
        return DataArray(Int8,nobs) #Vector{Union{Missing,Int8}}(nobs)
    end

    error(vtype, " is not a valid variable type in Stata.")
end

function alloc_unionvec(vtype,vfmt,nobs::Int64)

    # create a DataArray for the relevant type
    if 0 <= vtype < 2045 || vtype == 32768 # string variable
        return Vector{Union{Missing,String}}(nobs)
    elseif vtype == 65526
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return Vector{Union{Missing,DateTime}}(nobs)
        else
            return Vector{Union{Missing,Float64}}(nobs)
        end
    elseif vtype == 65527
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return Vector{Union{Missing,DateTime}}(nobs)
        else
            return Vector{Union{Missing,Float32}}(nobs)
        end
    elseif vtype == 65528
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return Vector{Union{Missing,DateTime}}(nobs)
        else
            return Vector{Union{Missing,Int32}}(nobs)
        end
    elseif vtype == 65529
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(nobs)
        else
            return Vector{Union{Missing,Int16}}(nobs)
        end
    elseif vtype == 65530
        return Vector{Union{Missing,Int8}}(nobs)
    end

    error(vtype, " is not a valid variable type in Stata.")
end

function memory_saved(da::DataArray)
    if sum(da.na) == length(da)
        return false
    end

    # total bytes
    tsize = sum(length.(dropna(s)))

    # reduced size
    lev = levels(dropna(s))
    if length(lev) > typmax(Int16)
        return false
    end
    rsize = sum(length.(lev)) + 4*length(da)

    if tsize > rsize
        return true
    end

    return false
end

function da2ca(da::DataArray)
    return CategoricalArray([isna(x) ? missing : x for x in da])
end

function read_stata(fn; categorize=true,verbose=false,exclude=[])
    df = DataFrame()
    label = Dict()

    read_stata!(fn, df, label, categorize=categorize,verbose=verbose,exclude=exclude)

    return (df,label)
end

function read_stata!(fn::String,df::DataFrame,label::Dict; categorize=true, verbose=false, exclude=[])

    fh = open(fn)

    # dta file (<stata_dta><header><release>)
    header = String(read(fh,UInt8,67))
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
    skip(fh,3)
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
        dslabel = String(read(fh,UInt8,dslabel_len))
    end

    # time stamp
    skip(fh,19)
    timestamplen = read(fh,Int8)
    if timestamplen > 0
        timestamp = String(read(fh,UInt8,timestamplen))
    end

    # map
    skip(fh,26) # </timestamp></header><map>
    statamap = Array{Int64}(14)
    for i in 1:14
        statamap[i] = read(fh,Int64)
    end

    # variable types
    skip(fh,22) # </map><variable_types>
    typelist = Array{UInt16}(nvar)
    for i in 1:nvar
        typelist[i] = read(fh,UInt16)
    end
    #print(Int.(typelist))

    # variable names
    skip(fh,27)
    varlist = Array{Symbol,1}(nvar)
    for i in 1:nvar
        varlist[i] = Symbol(strtonull(String(read(fh,UInt8,len_varname))))
    end

    # sort list
    skip(fh,21) # </varnames><sortlist>
    for i in 1:nvar
        srt = read(fh,Int16)
    end

    # formats
    skip(fh,22) # </sortlist><formats> + 2 (2 bytes left over from the previous sequence)
    fmtlist = Array{String,1}(nvar)
    for i in 1:nvar
        fmtlist[i] = strtonull(String(read(fh,UInt8,len_format)))
    end

    # value label names
    skip(fh,29) # </formats><value_label_names>
    labelnames = Array{String,1}(nvar)
    numvlabels = 0
    for i in 1:nvar
        labelnames[i] = strtonull(String(read(fh,UInt8,len_labelname)))

        # count the number of value labels
        if length(labelnames[i]) > 0
            numvlabels += 1
        end
    end

    # variable labels
    skip(fh,37) # </value_label_names><variable_labels>
    varlabels = Array{String}(nvar)
    for i in 1:nvar
        varlabels[i] = strtonull(String(read(fh,UInt8,len_varlabel)))
    end

    # characteristics - we will not import them
    # read until we hit '</characteristics>'
    skip(fh,35) # </variable_labels><characteristics>
    while (true)
        readuntil(fh,'<')
        if String(copy(read(fh,UInt8,5))) == "data>"
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

    # slurp the entire data section into memory
	io = IOBuffer()
    write(io,read(fh,UInt8,rlen*nobs))
	seek(io,0)

    # if there is strLs, read them now
    skip(fh,7) # </data>
    tst = String(read(fh,UInt8,7))
    if tst == "<strls>"
        # read strLs
        #
        # strL. Stata 13 introduced long strings up to 2 billon characters. strLs are
        # sperated by "GSO".
        # (v,o): Position in the data.frame.
        # t:     129/130 defines whether or not the strL is stored with a binary 0.
        # len:   length of the strL.
        # strl:  long string.
        strls = OrderedDict()
        t = OrderedDict()
        while (String(read(fh,UInt8,3)) == "GSO")
            v = read(fh,Int32)
            if release == 117
                o = read(fh,Int32)
            elseif release == 118
                o = read(fh,Int64)
            end
            t[(v,o)] = read(fh,UInt8)
            len = read(fh,UInt32)
            strls[(v,o)] = String(read(fh,UInt8,len))
        end
    end

    # read value labels
	#println(bytestring(readbytes(fh,41))) # </data><strls></strls><value_labels><lbl> ... </lbl>
    tst = String(read(fh,UInt8,5))
    if tst == "trls>"
        skip(fh,14)
    else
        error("Wrong position")
    end

    # loop through all value labels
    # define a Dict first
    value_labels = Dict()

    # numvlabels is defined above at the header section
    for i in 1:numvlabels
        #skip(fh,5) # <lbl>
		skipstr = String(read(fh,UInt8,5))
		if skipstr != "<lbl>"
			break
		end
        len = read(fh,Int32)
        labname = strtonull(String(read(fh,UInt8,len_labelname)))

        skip(fh,3) # padding
        numvalues = read(fh,Int32) # number of entries
        txtlen = read(fh,Int32) # length of value label text
        value_labels[labname] = OrderedDict()

        #
        offset = read(fh,Int32,numvalues) # offset
        values = read(fh,Int32,numvalues) # values
        valtext = String(read(fh,UInt8,txtlen)) # text table
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

        if verbose == true
            println("Processing variable ",j," ", varlist[j])
        end

        #df[varlist[j]] = alloc_dataarray(typelist[j],fmtlist[j],nobs)
        #df[varlist[j]] = alloc_unionvec(typelist[j],fmtlist[j],nobs)
        #da = alloc_dataarray(typelist[j],fmtlist[j],nobs)
        da = alloc_unionvec(typelist[j],fmtlist[j],nobs)

		for i in 1:nobs

            seek(io,numskip[j] + (i-1)*rlen)

            if 0 <= typelist[j] < 2045
                da[i] = strtonull(String(read(io,UInt8,typelist[j])))
                # if empty string, return NA
                if da[i] == ""
                    da[i] = missing
                end
            elseif typelist[j] == 32768 # long string
                if release == 117
                    v = read(io,Int32)
                    o = read(io,Int32)
                elseif release == 118
                    z = read(io,UInt8,8)
                    v = reinterpret(Int16,z[1:2])[1]
                    o = (reinterpret(Int64,z)[1] >> 16)
                end
                if (v,o) == (0,0)
                    da[i] = missing
                else
                    da[i] = strls[(v,o)]
                end
            elseif typelist[j] == 65526
                dataitemf64 = read(io,Float64)
                if dataitemf64 > 8.9884656743e307
                    da[i] = missing
                elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
                    # convert it to Julia date
                    da[i] = Date(1960,1,1) + Dates.Day(round(Int,dataitemf64))
                elseif fmtlist[j][1:3] == "%tc" || fmtlist[j][1:3] == "%tC"
                    da[i] = DateTime(1960,1,1,0,0,0) + Dates.Millisecond(round(Int,dataitemf64))
                else
                    da[i] = dataitemf64
                end
            elseif typelist[j] == 65527
                dataitemf32 = read(io,Float32)
                if dataitemf32 > 1.70141173319e38
                    da[i] = missing
                elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
                    # convert it to Julia date
                    da[i] = Date(1960,1,1) + Dates.Day(round(Int,dataitemf32))
                elseif fmtlist[j][1:3] == "%tc" || fmtlist[j][1:3] == "%tC"
                    da[i] = DateTime(1960,1,1,0,0,0) + Dates.Millisecond(round(Int,dataitemf32))
                else
                    da[i] = dataitemf32
                end
            elseif typelist[j] == 65528
                dataitemi32 = read(io,Int32)
                if dataitemi32 > 2147483620
                    da[i] = missing
                elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
                    # convert it to Julia date
                    da[i] = Date(1960,1,1) + Dates.Day(dataitemi32)
                elseif fmtlist[j][1:3] == "%tc" || fmtlist[j][1:3] == "%tC"
                    da[i] = DateTime(1960,1,1,0,0,0) + Dates.Millisecond(dataitemi32)
                else
                    da[i] = dataitemi32
                end
            elseif typelist[j] == 65529
                dataitemi16 = read(io,Int16)
                if dataitemi16 > 32740
                    da[i] = missing
                elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
                    # convert it to Julia date
                    da[i] = Date(1960,1,1) + Dates.Day(dataitemi16)
                else
                    da[i] = dataitemi16
                end
            elseif typelist[j] == 65530
                dataitemi8 = read(io,Int8)
                if dataitemi8 > 100
                    da[i] = missing
                else
                    da[i] = dataitemi8
                end
            end
        end

        # strls will be converted to categorical regardless of `categorize` option
        if typelist[j] == 32768
            #categorical!(df,varlist[j])
            df[varlist[j]] = categorial(da)
        # string variables
        # elseif 0 < typelist[j] < 2045
        #     # string variables can optionally be converted to categorical with the "categorize" option
        #     if categorize && in(varlist[j],exclude) && memory_saved(da)
        #         df[varlist[j]] = categorical(da)
        #     end
        # if an integer type and value labels exist, convert it to a categorical array and attach value labels
        elseif labelnames[j] != "" && in(typelist[j],[65528,65529,65530])
            value_labels[labelnames[j]][Missings.missing] = ""
            # df[varlist[j]] = CategoricalArrays.recode(da2ca(da),value_labels[labelnames[j]]...)
            df[varlist[j]] = CategoricalArrays.recode(categorical(da),value_labels[labelnames[j]]...)
        else
            df[varlist[j]] = da
        end
    end

    # variable labels
	variable_dict = Dict()

    # display formats
    format_dict = Dict()

    # value labels
    labname_dict = Dict()
    for i in 1:nvar
        variable_dict[varlist[i]] = varlabels[i]
        format_dict[varlist[i]] = fmtlist[i]
        labname_dict[varlist[i]] = labelnames[i]
    end

    label["variable"] = variable_dict
    label["format"] = format_dict
    label["label"] = labname_dict
    label["value"] = value_labels

	# close the file
	close(fh)

    # release memory
    # io = IOBuffer()
    # strls = Dict()
    # gc()

end
