#######################################################################
# Convert Stata DTA file to Julia DF
#######################################################################
"""
	read_stata(fn::String; chunks::Int=10)

converts a stata datafile `fn` to Julia DataFrame. An original data file bigger than 100MB will be read in `chunks` (default = 10)
to save memory. 
"""
function read_stata(fn::String; chunks::Int=10)

    fh = open(fn, "r")

    # dta file (<stata_dta><header><release>)
    header = String(read(fh, 67))
    if header[2:10] != "stata_dta"
        error("Not a version 13 or 14 data file")
    end

    # data format version
    release = parse(Int16, header[29:31])
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
        error("Can't convert data format version ", release, ".")
    end

    # byte order: LSF or MSF
    byteorder = header[53:55]
    if byteorder == "MSF"
        error("Big-endian data are not supported.")
    end

    # number of variables
    skip(fh, 3) # <K>
    nvar = read(fh, Int16)

    # number of observations
    skip(fh, 7) #</K><N>
    nobs = Int(read(fh, Int32))

    # dataset label length
    if release == 117
        skip(fh, 11)
        dslabel_len = read(fh, Int8)
    elseif release == 118
        skip(fh, 15)
        dslabel_len = read(fh, Int16)
    end

    # read the label
    dslabel = ""
    if dslabel_len > 0
        dslabel = String(read(fh, dslabel_len))
    end

    # time stamp
    skip(fh, 19)
    timestamplen = read(fh, Int8)
    if timestamplen > 0
        timestamp = String(read(fh, timestamplen))
    end

    # map
    skip(fh, 26) # </timestamp></header><map>
    statamap = Vector{Int64}(undef, 14)
    read!(fh, statamap)

    # variable types
    skip(fh, 22) # </map><variable_types>
    typelist = Vector{UInt16}(undef, nvar)
    read!(fh, typelist)

    # variable names
    skip(fh, 27)
    varlist = Vector{Symbol}(undef, nvar)
    for i in 1:nvar
        varlist[i] = Symbol(strtonull(String(read(fh, len_varname))))
    end

    # sort list
    skip(fh, 21) # </varnames><sortlist>
    for i in 1:nvar
        srt = read(fh, Int16)
    end

    # formats
    skip(fh, 22) # </sortlist><formats> + 2 (2 bytes left over from the previous sequence)
    fmtlist = Vector{String}(undef, nvar)
    for i in 1:nvar
        fmtlist[i] = strtonull(String(read(fh, len_format)))
    end

    # value label names
    skip(fh, 29) # </formats><value_label_names>
    valuelabels = Vector{String}(undef, nvar)
    numvlabels = 0
    for i in 1:nvar
        valuelabels[i] = strtonull(String(read(fh, len_labelname)))

        # count the number of value labels
        if length(valuelabels[i]) > 0
            numvlabels += 1
        end
    end

    # variable labels
    skip(fh, 37) # </value_label_names><variable_labels>
    varlabels = Vector{String}(undef, nvar)
    for i in 1:nvar
        varlabels[i] = strtonull(String(read(fh, len_varlabel)))
    end

    # characteristics - we will not import them
    # read until we hit '</characteristics>'
    skip(fh, 35) # </variable_labels><characteristics>
    while (true)
        readuntil(fh, '<')
        if String(copy(read(fh, 5))) == "data>"
            break
        end
    end

    # nubmer of bytes for each variable
    numbytes = get_numbytes(typelist, nvar)

    # total length of each observation
    rlen = sum(numbytes)

    # number of bytes to skip in IOBuffer
    numskip = zeros(Int, nvar)
    numskip[1] = 0
    for i in 2:length(numbytes)
        numskip[i] = numbytes[i-1] + numskip[i-1]
    end

    # save the start position of the data section
    data_pos = position(fh)

    # skip the data section for now
    skip(fh, rlen * nobs)

    # if there is strLs, read them now
    skip(fh, 7) # </data>
    tst = String(read(fh, 7))
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
        while (String(read(fh, 3)) == "GSO")
            v = read(fh, Int32)
            if release == 117
                o = read(fh, Int32)
            elseif release == 118
                o = read(fh, Int64)
            end
            t[(v, o)] = read(fh, UInt8)
            len = read(fh, UInt32)
            strls[(v, o)] = String(read(fh, len))
        end
    end

    tst = String(read(fh, 5))
    if tst == "trls>"
        skip(fh, 14)
    else
        error("Wrong position")
    end

    # read value labels
    # loop through all value labels
    # define a Dict first

    value_labels = Dict()

    # numvlabels is defined above at the header section
    for i in 1:numvlabels

        skipstr = String(read(fh, 5))
        if skipstr != "<lbl>"
            break
        end
        len = read(fh, Int32)
        labname = Symbol(strtonull(String(read(fh, len_labelname))))

        skip(fh, 3) # padding
        numvalues = read(fh, Int32) # number of entries
        txtlen = read(fh, Int32) # length of value label text
        value_labels[labname] = Dict()

        #
        offset = Vector{Int32}(undef, numvalues)
        read!(fh, offset) # offset
        values = Vector{Int32}(undef, numvalues)
        read!(fh, values) # values
        valtext = String(read(fh, txtlen)) # text table

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
        skip(fh, 6) # </lbl>
    end

    variable_dict = Dict()
    lblname_dict = Dict()
    for i in 1:nvar
        if length(varlabels[i]) > 0
            variable_dict[varlist[i]] = varlabels[i]
        end

        if haskey(value_labels, i)
            lblname_dict[varlist[i]] = Symbol(value_labels[i])
        end
    end

    # read data now
    seek(fh, data_pos)

    # if the data size < 100MB, then
    # slurp the entire data section into memory
    # otherwise, we will read the data by smaller batches
    io = IOBuffer()
    if rlen * nobs < 100_000_000
        write(io, read(fh, rlen * nobs))
        seek(io, 0)
        rdf = _read_dta(io, release, rlen, nobs, nvar, varlist, variable_dict, typelist, fmtlist, lblname_dict, value_labels, numskip, strls)
    else
        len = max(100000, ceil(Int, nobs / chunks))
        totlen = nobs
        rdf = DataFrame()
        for nread in 1:ceil(Int, nobs / len)
            if len < totlen
                totlen -= len
            else
                len = totlen
            end
            seek(io, 0)
            write(io, read(fh, rlen * len))
            seek(io, 0)

            rdf = vcat(rdf, _read_dta(io, release, rlen, len, nvar, varlist, variable_dict, typelist, fmtlist, lblname_dict, value_labels, numskip, strls))
        end
    end

    # close the file
    close(fh)

    return rdf
end

# function read_labels(fn::String)
# 	return read_stata(fn,read_labels=true)
# end
function _read_dta(io, release, rlen, len, nvar, varlist, varlabels, typelist, fmtlist, lblname, vallabels, numskip, strls)

    df = DataFrame()

    dataitemf32::Float32 = 0.0
    dataitemf64::Float64 = 0.0
    dataitemi8::Int8 = 0
    dataitemi16::Int16 = 0
    dataitemi32::Int32 = 0
    v::Int32 = 0
    o::Int64 = 0
    z = zeros(UInt8, 8)

    # interate over the number of variables
    for j in 1:nvar

        df[!, varlist[j]] = alloc_array(typelist[j], fmtlist[j], len)

        # len == number of observation in the batch
        for i in 1:len

            seek(io, numskip[j] + (i - 1) * rlen)

            if 0 <= typelist[j] < 2045
                df[i, j] = strtonull(String(read(io, typelist[j])))
                # if empty string, return missing
                if df[i, j] == ""
                    df[i, j] = missing
                end
            elseif typelist[j] == 32768 # long string
                if release == 117
                    v = read(io, Int32)
                    o = read(io, Int32)
                elseif release == 118
                    z = read(io, 8)
                    v = reinterpret(Int16, z[1:2])[1]
                    o = (reinterpret(Int64, z)[1] >> 16)
                end
                if (v, o) == (0, 0)
                    df[i, j] = missing
                else
                    df[i, j] = strtonull(strls[(v, o)])
                end
            elseif typelist[j] == 65526
                dataitemf64 = read(io, Float64)
                if dataitemf64 > 8.9884656743e307
                    df[i, j] = missing
                elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
                    # convert it to Julia date
                    df[i, j] = Date(1960, 1, 1) + Dates.Day(round(Int, dataitemf64))
                elseif fmtlist[j][1:3] == "%tc" || fmtlist[j][1:3] == "%tC"
                    df[i, j] = DateTime(1960, 1, 1, 0, 0, 0) + Dates.Millisecond(round(Int, dataitemf64))
                else
                    df[i, j] = dataitemf64
                end
            elseif typelist[j] == 65527
                dataitemf32 = read(io, Float32)
                if dataitemf32 > 1.70141173319e38
                    df[i, j] = missing
                elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
                    # convert it to Julia date
                    df[i, j] = Date(1960, 1, 1) + Dates.Day(round(Int, dataitemf32))
                elseif fmtlist[j][1:3] == "%tc" || fmtlist[j][1:3] == "%tC"
                    df[i, j] = DateTime(1960, 1, 1, 0, 0, 0) + Dates.Millisecond(round(Int, dataitemf32))
                else
                    df[i, j] = dataitemf32
                end
            elseif typelist[j] == 65528
                dataitemi32 = read(io, Int32)
                if dataitemi32 > 2147483620
                    df[i, j] = missing
                elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
                    # convert it to Julia date
                    df[i, j] = Date(1960, 1, 1) + Dates.Day(dataitemi32)
                elseif fmtlist[j][1:3] == "%tc" || fmtlist[j][1:3] == "%tC"
                    df[i, j] = DateTime(1960, 1, 1, 0, 0, 0) + Dates.Millisecond(dataitemi32)
                else
                    df[i, j] = dataitemi32
                end
            elseif typelist[j] == 65529
                dataitemi16 = read(io, Int16)
                if dataitemi16 > 32740
                    df[i, j] = missing
                elseif fmtlist[j] == "%d" || fmtlist[j][1:3] == "%td"
                    # convert it to Julia date
                    df[i, j] = Date(1960, 1, 1) + Dates.Day(dataitemi16)
                else
                    df[i, j] = dataitemi16
                end
            elseif typelist[j] == 65530
                dataitemi8 = read(io, Int8)
                if dataitemi8 > 100
                    df[i, j] = missing
                else
                    df[i, j] = dataitemi8
                end
            end
        end
        # strls will be converted to categorical regardless of `categorize` option
        if typelist[j] == 32768
            categorical!(df, varlist[j])
        end

        # for integer variables that have formats
        # convert them into CategoricalArrays with the appropriate value labels
        if typelist[j] in (65528, 65529, 65530) && haskey(lblname, varlist[j])
            values!(df, varlist[j], vallabels[lblname[varlist[j]]])
        end

        # variable label
        if haskey(varlabels, varlist[j])
            TableMetadataTools.label!(df, varlist[j], varlabels[varlist[j]])
        end

        # for vectors without missing values, convert the vector to an appropirate type
        if sum(ismissing.(df[!, varlist[j]])) == 0
            df[!, varlist[j]] = convert(Vector{nonmissingtype(eltype(df[!, varlist[j]]))}, df[!, varlist[j]])
        end
    end

    return Stella.dfcompress(df)
end

function get_numbytes(typelist, nvar)
    nb = Vector{UInt16}(undef, nvar)
    for i in 1:nvar
        if 0 < typelist[i] < 2045
            nb[i] = typelist[i]
        elseif typelist[i] == 32768
            nb[i] = 8
        elseif typelist[i] == 65526
            nb[i] = 8
        elseif typelist[i] == 65527
            nb[i] = 4
        elseif typelist[i] == 65528
            nb[i] = 4
        elseif typelist[i] == 65529
            nb[i] = 2
        elseif typelist[i] == 65530
            nb[i] = 1
        end
    end
    return nb
end

function strtonull(str::String)
    n = findfirst('\0', str)
    n == nothing && return str
    return str[1:n-1]
end

function alloc_array(vtype, vfmt, nobs::Int64)

    if 0 <= vtype < 2045 || vtype == 32768 # string variable
        return Vector{Union{Missing,String}}(missing, nobs)
    elseif vtype == 65526
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(missing, nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return Vector{Union{Missing,DateTime}}(missing, nobs)
        else
            return Vector{Union{Missing,Float64}}(missing, nobs)
        end
    elseif vtype == 65527
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(missing, nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return Vector{Union{Missing,DateTime}}(missing, nobs)
        else
            return Vector{Union{Missing,Float32}}(missing, nobs)
        end
    elseif vtype == 65528
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(missing, nobs)
        elseif vfmt[1:3] == "%tc" || vfmt[1:3] == "%tC"
            return Vector{Union{Missing,DateTime}}(missing, nobs)
        else
            return Vector{Union{Missing,Int32}}(missing, nobs)
        end
    elseif vtype == 65529
        if vfmt == "%d" || vfmt[1:3] == "%td"
            return Vector{Union{Missing,Date}}(missing, nobs)
        else
            return Vector{Union{Missing,Int16}}(missing, nobs)
        end
    elseif vtype == 65530
        return Vector{Union{Missing,Int8}}(missing, nobs)
    end

    error(vtype, " is not a valid variable type in Stata.")
end


#######################################################################
# Convert Julia DF to Stata DTA format
#######################################################################

missingval = Dict(
    65530 => 101,
    65529 => 32_741,
    65528 => 2_147_483_621,
    65527 => typemax(Float32),
    65526 => typemax(Float64)
)

vtype = Dict(
    Bool => 65530,
    Int8 => 65530,
    Int16 => 65529,
    Int32 => 65528,
    Float32 => 65527,
    Float64 => 65526,
    Date => 65528,
    DateTime => 65526
)

function write_stata(fn::String,outdf::AbstractDataFrame; maxbuffer = 10_000_000, verbose = false)

    # open the output dataset file
    if fn[end-3:end] != ".dta"
        fn = string(fn, ".dta")
    end

    outdta = open(fn,"w")

    # allowed variable types
    notallowed = [ in(x, [Bool, Int8, Int16, Int32, Int64, Float32, Float64, String, Date, DateTime]) ? false : true for x in dtypes(outdf)]

    # empty variables
    allmiss = [ sum(ismissing.(x)) == size(outdf,1) ? true : false for x in eachcol(outdf)]

    # large Int64 values that cannot be saved as Int32
    lint64 = falses(ncol(outdf))
    for i in 1:ncol(outdf)
        if eltype2(outdf[:,i]) == Int64 && allmiss[i] == false
            tvec = collect(skipmissing(outdf[:,i]))
            if maximum(tvec) > 2_147_483_620 || minimum(tvec) < −2_147_483_647
                lint64[i] = true
            end
        end
    end
    
    # subset
    df = outdf[:,findall(x->x == true, [ notallowed[x] || allmiss[x] || lint64[x] ? false : true for x in 1:ncol(outdf)])]

    # report exclusions
    if verbose
        println("These variables will NOT be exported because of their data types:\n")
        for (i, v) in enumerate(names(outdf))
            notallowed[i] && println(@sprintf("%-30s\t%-20s",v, eltype2(outdf[:,v])))
        end
        println("\n\nThese variables will be excluded variables because they are empty (100% missing).\n")
        for (i, v) in enumerate(names(outdf))
            allmiss[i] && println(@sprintf("%-30s\t%-20s",v, eltype2(outdf[:,v])))
        end
        println("\n\nThese variables will be excluded variables because they are Int64 variables that contain values larger than Int32 can hold.\n")
        for (i, v) in enumerate(names(outdf))
            lint64[i] && println(@sprintf("%-30s\t%-20s",v, eltype2(outdf[:,v])))
        end
    end

    # data types
    datatypes = dtypes(df)

    # cols and rows
    (rows, cols) = size(df)

    # -----------------------------------------------------
    # release 118 format parameters
    len_varname = 129
    len_format = 57
    len_labelname = 129
    len_varlabel = 321

    # header
    write(outdta,"<stata_dta><header><release>118</release><byteorder>LSF</byteorder><K>")

    # number of variables
    write(outdta,Int16(cols))
    write(outdta,"</K><N>")

    # number of observations
    write(outdta, Int64(rows))
    write(outdta, "</N><label>")

    # assume no data label
    write(outdta,UInt16(0))
    write(outdta,"</label><timestamp>")
    
    # timestamp
    ts = Dates.format(now(), "dd uuu yyyy HH:MM")
    write(outdta,UInt8(length(ts)))
    write(outdta,string(ts,"</timestamp></header>"))

    # -----------------------------------------------------
    # map
    m = zeros(Int64, 14)
    m[2] = Int64(position(outdta))
    write(outdta,"<map>")
    write(outdta,m)
    write(outdta,"</map>")

    # -----------------------------------------------------
    # variable types
    (typelist, numbytes) = get_types(df)
    m[3] = Int64(position(outdta))
    write(outdta,"<variable_types>")
    write(outdta,UInt16.(typelist))
    write(outdta,"</variable_types>")

    # variable names
    m[4] = Int64(position(outdta))
    write(outdta,"<varnames>")
    write(outdta, get_varnames(df,len_varname))
    write(outdta,"</varnames>")

    # sortlist
    m[5] = Int64(position(outdta))
    write(outdta,"<sortlist>")
    write(outdta,zeros(Int16, cols+1))
    write(outdta,"</sortlist>")

    # formats - no formats for now
    m[6] = Int64(position(outdta))
    write(outdta,"<formats>")
    write(outdta,get_formats(df,typelist,len_format))
    write(outdta,"</formats>")

    # value label names
    m[7] = Int64(position(outdta))
    write(outdta,"<value_label_names>")
    write(outdta,get_label_names(df,len_labelname))
    write(outdta,"</value_label_names>")

    # variable labels
    m[8] = Int64(position(outdta))
    write(outdta,"<variable_labels>")
    write(outdta,get_variable_labels(df,len_varlabel))
    write(outdta,"</variable_labels>")

    # characteristics - nothing to output
    m[9] = Int64(position(outdta))
    write(outdta,"<characteristics></characteristics>")

    # total length of a row
    rlen = sum(numbytes)
    
    # ---------------------------------------------------------
    # the rest of the map section data
    m[10] = Int64(position(outdta))
    write(outdta,"<data>")

    # --------------------------------------------------------=
    # combine rows into one iobuffer and write
    chunks = ceil(Int32, rlen * rows / maxbuffer)
    nobschunk = chunks == 1 ? nobschunk = rows : ceil(Int32, rows / (chunks - 1))
    for i = 1:chunks
        from = 1 + (i-1)*nobschunk
        to = min(from + nobschunk - 1, rows)
        write(outdta,write_chunks(@view(df[from:to,:]), datatypes, typelist, rlen))
    end
    write(outdta,"</data>")

    # strL section - no strLs
    m[11] = Int64(position(outdta))
    write(outdta, "<strls></strls>")

    # value labels
    m[12] = Int64(position(outdta))
    write(outdta,"<value_labels>")
    write(outdta, get_value_labels(df))
    write(outdta,"</value_labels>")

    # at the end of stata_dta section
    m[13] =  Int64(position(outdta))

    # end of file
    write(outdta,"</stata_dta>")
    m[14] = Int64(position(outdta)) 

    # seek back to the map section and overwrite the map data
    seek(outdta, m[2]+length("<map>"))
    write(outdta, m)

    # flush iostream and close
    flush(outdta)
    close(outdta)

end

function write_chunks(outdf, datatypes, typelist, rlen)
    iobuf = IOBuffer()
    for (k,dfrow) in enumerate(eachrow(outdf))
        loc = position(iobuf)
        for (i,v) in enumerate(dfrow)
            if isa(outdf[:,i], CategoricalArray) # CategoricalArray
                if eltype2(outdf[:,i]) == String
                    write(iobuf, Int32(ismissing(v) ? missingval[Int32] : outdf[:,i].pool.invindex[v]))
                # else
                #     write(iobuf, datatypes[i](ismissing(v) ? missingval[typelist[i]] : unwrap(v)))
                end
            elseif datatypes[i] == String
                write(iobuf, ismissing(v) ? repeat('\0', typelist[i]) : string(v, repeat('\0', typelist[i] - sizeof(v))))
            elseif datatypes[i] == Date
                write(iobuf, Int32(ismissing(v) ? missingval[65528] : Dates.value(v - Date(1960,1,1))))
            elseif datatypes[i] == DateTime
                write(iobuf, Float64(ismissing(v) ? missingval[65526] : Dates.value(v - DateTime(1960,1,1))))
            elseif datatypes[i] == Bool
                write(iobuf, Int8(ismissing(v) ? missingval[65530] : v == true ? 1 : 0))
            elseif datatypes[i] == Int64
                write(iobuf, Int32(ismissing(v) ? missingval[65528] : v))
            else
                write(iobuf, datatypes[i](ismissing(v) ? missingval[typelist[i]] : v))
            end
        end
        if position(iobuf) - loc != rlen
            error("Data overrun on observation ",k, "; record length = ", rlen,"; current length = ", position(iobuf) - loc)
        end

    end
    return take!(iobuf)
end

function dtypes(outdf)
    t = []
    for i in 1:ncol(outdf)
        if isa(outdf[:,i], CategoricalArray)
            typ = eltype2(outdf[:,i])
            if typ == String
                push!(t, Int32)
            else
                push!(t, typ)
            end
        else
            push!(t, eltype2(outdf[:,i]))
        end
    end
    return t
end

function getmaxbytes(s::AbstractArray)
    if isa(s, CategoricalArray) && eltype2(s) <: CategoricalString
        return maximum(sizeof.(levels(s)))
    end

    if nmissing(s) == size(s, 1)
        return 0
    end

    return maximum(sizeof.(skipmissing(s)))
end

function get_types(outdf)

    bytesize = Dict(
        65526 => 8,
        65527 => 4,
        65528 => 4,
        65529 => 2,
        65530 => 1,
    )

    tlist = zeros(Int32,ncol(outdf))
    numbytes = zeros(Int32,ncol(outdf))
    for i in 1:ncol(outdf)
        if isa(outdf[:,i], CategoricalArray)
            typ = eltype2(outdf[:,i])
            if typ == String
                tlist[i] = 65528
                numbytes[i] = 4
            elseif typ == Bool
                tlist[i] = 65530
                numbytes[i] = 1
            else
                tlist[i] = vtype[typ]
                numbytes[i] = bytesize[tlist[i]]
            end
        else
            typ = eltype2(outdf[:,i])
            if haskey(vtype,typ)
                tlist[i] = vtype[typ]
                numbytes[i] = bytesize[tlist[i]]
            elseif typ == Bool
                tlist[i] = 65530
                numbytes[i] = 1
            elseif typ == Int64
                tlist[i] = vtype[Int32]
                numbytes[i] = bytesize[tlist[i]]
            elseif typ == String
                maxlen = Stella.getmaxbytes(outdf[:,i])
                if maxlen < 2045
                    tlist[i] = maxlen
                    numbytes[i] = tlist[i]
                else
                    tlist[i] = 32768
                    numbytes[i] = 4
                end
            end
        end
    end

    return (tlist, numbytes)
end

function get_varnames(outdta,len)

    varstring = String[]
    for v in names(outdta)
        vlen = sizeof(v)
        if vlen < len - 1
            push!(varstring, string(v,repeat('\0',len - vlen)))
        else
            push!(varstring, string(v[1:end-1],'\0'))
        end
    end
    return join(varstring,"")
end

function get_formats(outdf,typelist,len)

    fvec = String[]
    for i in 1:ncol(outdf)
        if typelist[i] < 2045
            fmt = string("%-",typelist[i],"s")
            push!(fvec,string(fmt, repeat('\0',len - sizeof(fmt))))
        elseif typelist[i] == 65528 && eltype2(outdf[:,i]) == Date
            push!(fvec,string("%tdNN-DD-CCYY",repeat('\0',len - 13)))
        elseif typelist[i] in (65528,65529,65530)
            push!(fvec,string("%8.0g",repeat('\0',len - 5)))
        elseif typelist[i] == 65527 # float
            push!(fvec,string("%6.2f",repeat('\0',len - 5)))
        elseif typelist[i] == 65526 && eltype2(outdf[:,i]) == DateTime
            push!(fvec,string("%tc",repeat('\0',len - 3)))
        elseif typelist[i] == 65526
            push!(fvec,string("%11.1f",repeat('\0',len - 6)))
        end
    end    

    return join(fvec,"")
end

function get_label_names(outdf,len)

    lvec = String[]
    for i in 1:size(outdf,2)
        if isa(outdf[:,i], CategoricalArray) && eltype(levels(outdf[:,i])) == String
            lblname = string("fmt",i)
            push!(lvec,string(lblname, repeat('\0',len - sizeof(lblname))))
        else
            push!(lvec,repeat('\0',len))
        end
    end    
    return join(lvec,"")
end

function get_variable_labels(outdf, len)
    lbls = labels(outdf)
    for i in 1:length(lbls)
        lbls[i] = string(lbls[i], repeat('\0', len - sizeof(lbls[i])))
    end
    return join(lbls,"")
end

function get_value_labels(outdf)

    iobf = IOBuffer()
    for (j,v) in enumerate(propertynames(outdf))
        if isa(outdf[:,v], CategoricalArray) && eltype2(outdf[:,v]) == String

            invindex = outdf[:,v].pool.invindex
            n = length(invindex)
            off = zeros(Int32,n)
            val = Int32.(values(invindex))
            txt = ""
            for (i,vv) in enumerate(keys(invindex))
                off[i] = sizeof(txt)
                txt = string(txt, vv, '\0')
            end
            txtlen = sizeof(txt)
            len = 8 + 8*n + txtlen

            # write the value label to the iobuffer
            write(iobf,"<lbl>")
            write(iobf,Int32(len))
            fmt = string("fmt",j)
            write(iobf,string(fmt,repeat('\0',129 - sizeof(fmt)), "   ")) # 3 spaces padded
            write(iobf,Int32(n))
            write(iobf,Int32(txtlen))
            write(iobf,off)
            write(iobf,val)
            write(iobf,txt)
            write(iobf,"</lbl>")
        end
    end

    return take!(iobf)
end

