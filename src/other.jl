"""
    identify_condition(df::AbstractDataFrame,cols::Vector{Symbol},codes::Vector{String})

produces a vector of Boolean `true` or `false` depending on whether columns `cols` contain
any of the string values in `codes`. This function is designed to easily identify comorbid
conditions in inpatient or outpatient records that are coded in ICD-9 or ICD-10 codes.
The values in `cols` are compared to the values (`code`) in `codes` only up to the length of
the `code`. For example, when detecting "diabetes", one may use "250.xx" in ICD-9,
which can be specified as identify_condition(df,[:diag1,:diag2,:diag3,:diag4,:diag5],["250"]).
All values that start with "250" in :diag1 - :diag4 (e.g., "25010") will match.
"""
function identify_condition(ip::AbstractDataFrame,vicd::Vector{Symbol},codes::Vector{String})

    retvec = zeros(Bool,nrow(ip))

    for i in 1:nrow(ip)
        for v in vicd
            if ismissing(ip[i,v])
                continue
            end

            vlen=length(ip[i,v])
            for c in codes
                clen = length(c)
                if vlen >= clen && c == ip[i,v][1:length(c)]
                    retvec[i] = true
                    break
                end
            end
            if retvec[i] == true
                break
            end
        end
    end
    return retvec
end

"""
    identify_condition2(df::AbstractDataFrame,cols::Vector{Symbol},codes::Vector{String})

produces a vector of Boolean `true` or `false` depending on whether columns `cols` contain
any of the string values in `codes`. This function is designed to easily identify comorbid
conditions in inpatient or outpatient records that are coded in ICD-9 or ICD-10 codes.
The values in `cols` are compared to the values in `codes` exactly as they are specified.
"""
function identify_condition2(ip::AbstractDataFrame,vicd::Vector{Symbol},codes::Vector{String})

    retvec = zeros(Bool,nrow(ip))

    for i in 1:nrow(ip)
        for v in vicd
            if ismissing(ip[i,v])
                continue
            end

            if ip[i,v] in codes
                retvec[i] = true
            end
        end
    end
    return retvec
end

using Survival

"""
    kaplanmeier(df, event, by = nothing)

Plots Kaplan-Meier estimates.
"""
function kaplanmeier(df, event, by=nothing)

    plt = nothing
    if by == nothing
        km = fit(KaplanMeier, df[:, event])
        plot(vcat(0, km.events.time), vcat(1.0, km.survival), linetype=:steppost, ylims=(0, 1))
        # return nothing
    else

        kvec = []

        for (i, v) in enumerate(sort(unique(skipmissing(df[:, by]))))
            df2 = filter(x -> !ismissing(x[by]) && x[by] == v, df)
            push!(kvec, fit(KaplanMeier, df2[!, event]))
            if i == 1
                plt = Plots.plot(vcat(0, kvec[1].events.time),
                    vcat(1.0, kvec[1].survival),
                    linetype=:steppost,
                    ylims=(0, 1),
                    xlabel="Analysis time",
                    ylabel="Survival estimates",
                    label=string(v))

            else
                Plots.plot!(plt, vcat(0, kvec[i].events.time),
                    vcat(1.0, kvec[i].survival),
                    linetype=:steppost,
                    ylims=(0, 1),
                    label=string(v))

            end
        end
    end
    return plt
end

function nelsonaalen(df, event, by=nothing)

    plt = nothing
    if by == nothing
        na = fit(NelsonAalen, df[:, event])
        plot(vcat(0, na.events.time), vcat(0, na.chaz), linetype=:steppost, ylims=(0, 1))
    else
        navec = []

        for (i, v) in enumerate(sort(unique(skipmissing(df[:, by]))))
            df2 = filter(x -> !ismissing(x[by]) && x[by] == v, df)
            push!(navec, fit(NelsonAalen, df2[!, event]))
            if i == 1
                plt = Plots.plot(vcat(0, navec[1].events.time),
                    vcat(0, navec[1].chaz),
                    linetype=:steppost,
                    ylims=(0, 1),
                    xlabel="Analysis time",
                    ylabel="Cumulative hazard rate",
                    label=string(v))
            else
                Plots.plot!(plt, vcat(0, navec[i].events.time),
                    vcat(0, navec[i].chaz),
                    linetype=:steppost,
                    ylims=(0, 1),
                    label=string(v))

            end
        end
    end
    return plt
end

# struct LogRank
#     observed::Vector{Int64}()
#     expected::Vector{Float64}()
#     nobs::Int64
#     dof::Int64
#     chi2::Float64
#     pvalue::Float64
# end

# function Base.show(io::IO, val::LogRank)
#     show(io, ev.time)
#     iscensored(ev) && print(io, '+')
#     return nothing
# end

"""
    logrank(df, event, by)

Performs log-rank test for the groups `by`.
"""
function logrank(df, event, by)

    # number of groups
    # ba = completecases(df[!,[event,by]])
    # df2 = df[ba,:]
    groups = sort(unique(df[!,by]))
    n_groups = length(groups)

    # perform Kaplan-Meier analysis
    km = Vector{KaplanMeier{Float64, Int64}}(undef,n_groups)
    times = zeros(Int64,n_groups)

    for (i,v) in enumerate(groups)
        km[i] = fit(KaplanMeier, df[ df[!,by] .== v, event])
        times[i] = km[i].events.time[end]
    end

    # lengths
    ntimes = maximum(times)

    # events
    events = zeros(Int64, n_groups, ntimes)
    for i in 1:n_groups
        events[i, km[i].events.time] .= km[i].events.nevents
    end

    # N at risk
    atrisk = zeros(Int64, n_groups, ntimes)
    for i in 1:n_groups
        atrisk[i, 1] = km[i].events.natrisk[1]
        
        for j in 2:ntimes
            atrisk[i, j] = atrisk[i, j-1] - events[i, j-1]
        end
    end

    # Observed events
    O = vec(sum(events, dims=1))

    # Total N at risk
    N = vec(sum(atrisk, dims=1))

    # Observed Rate
    Or = O ./ N

    # Expected values
    E = zeros(Float64, n_groups, ntimes)
    for i in 1:n_groups
        E[i, :] = Or .* vec(atrisk[i, :])
    end

    o = vec(sum(events, dims=2))
    e = vec(sum(E, dims=2))

    chi2 = sum((o .- e) .^ 2 ./ e)

    dof = (length(o) - 1) * (length(e) - 1)
    pval = ccdf(Chisq(dof), chi2)

    return (o, e, size(df,1), chi2, dof, pval)
end


"""
    st2ncc(::AbstractDataFrame, ev::EventTime; ncontrol = 1, matchvars=nothing)

Creates a nested case-control design from a survival data. It requires an EventTime
variable as defined in the Survival.jl package. The output data will be sorted
by :_set and :_case. Options are:

ncontrol - number of controls to find for each case
matchvars - a vector of matching factors

This function generates the following variables:

_case - 1 for case and 0 for controls
_set  - group (or set) ID
_nset - number of rows in the group (or set)
_time - time to event
"""
function st2ncc(df::AbstractDataFrame, ev; ncontrol=1, matchvars=nothing)

    # time
    df._time = [x[ev].time for x in eachrow(df)]

    # non-matches
    nonmatch = 0

    # empty output data
    dfout = DataFrame()

    # events
    dfev = df[isevent.(df[:, ev]), :]
    dfev._case = ones(Int8, nrow(dfev))
    dfev._set = collect(1:nrow(dfev))

    # iterate over all event rows
    sort!(dfev,[ev])
    for i = 1:nrow(dfev)

        # define a risk set consisting of those who have not experienced an event until the ev's eventtime
        df2 = filter(x -> x[ev] > dfev[i, ev], df)

        # if there are matchvars
        if matchvars != nothing
            for m in matchvars
                df2 = filter(x -> x[m] == dfev[i, m], df2)
            end
        end

        # no matches found
        if nrow(df2) == 0
            nonmatch += 1
            continue
        end

        # if df2 has more than ncontrol rows, select matches randomly
        if nrow(df2) > ncontrol
            df2 = df2[rand(1:nrow(df2), ncontrol), :]
        end

        df2._set .= dfev._set[i]
        df2._case .= 0
        df2._time .= dfev._time[i]
        dfout = vcat(dfout, df2)
    end

    # drop sets without any controls
    if nonmatch > 0
        println(nonmatch, " cases were dropped. No any matches were found.")
    end
    df2 = vcat(dfev, dfout)
    df3 = select(combine(groupby(df2, :_set), nrow => :_nset), [:_set, :_nset])
    df2 = leftjoin(df3[df3._nset .> 1,:], df2, on=:_set)

    # clean up and return
    return sort(df2, [:_set, :_case])
end

"""
    elixhauser9!(df, icdvars::Vector; drg = nothing)

Creates 30 variables in the `df` that indicates presence (1 or 0 otherwise) of 30 comorbidities in the
Elixhauser comorbidity index. This program is based on the SAS version of the AHRQ HCUP comorbidity
software version 3.6 using ICD-9 DX codes and MS-DRG V28. The program takes the following arguments:

- df - DataFrame that contains claims records

- icdvars - a vectors of variables of type Symbol or String that contain ICD-9 diagnostic codes. If you are computing comorbidities based on an inpatient dataset, you may want to omit the "principal diagnosis" variable and specify the `drg` variable (see below).

- drg - a variable name of type Symbol or String that contain DRG codes
"""
function elixhauser9!(df, icdvars::Vector; drg=nothing)

    # load data
    elixdata = JLD2.load(joinpath(@__DIR__, "..", "data", "elixhauser_v9.jld2"))
    icd9map = elixdata["icd9map"]
    condmap = elixdata["condmap"]
    drgmaptmp = elixdata["drgmap"]
    tmpmap = elixdata["tmpmap"]
    description = elixdata["description"]

    # drg - change datatype for DRG values according to the type in the DF
    if drg != nothing
        drgtype = nonmissingtype(eltype(df[:, drg]))
        if drgtype == String
            drgmap = drgmaptmp
        else
            drgmap = Dict(parse.(Int16, keys(drgmaptmp)) .=> values(drgmaptmp))
        end
    end

    # create comorbidity variables
    for i in 1:30
        df[:, condmap[i]] = zeros(Int8, nrow(df))
        label!(df, condmap[i], description[i])
    end

    for i in 1:nrow(df)

        # initialize tmpflags and htn variables
        tmpflg = falses(10)
        htn = 0
        htncx = 0

        # ICD-9 codes
        for icdvar in icdvars

            icd = df[i, icdvar]
            if ismissing(icd) || icd in ("", " ")
                continue
            end
            if haskey(icd9map, icd)
                # find the index for the ICD-9 code
                idx = icd9map[icd] # an ICD code can be mapped to 2 conditions
                vv = condmap[idx]
                df[i, vv] = 1
            end

            # temporary formats for HTNCX, CHF, RENLFAIL 
            if haskey(tmpmap, icd)
                tmpflg[tmpmap[icd]] = 1
            end
        end

        # temporary formats for HTNCX, CHF, RENLFAIL
        if any(tmpflg)
            htncx = 1
            if tmpflg[3] || tmpflg[7]
                df[i, :chf] = 1
            elseif tmpflg[5] || tmpflg[8]
                df[i, :renlfail] = 1
            elseif tmpflg[9]
                df[i, :chf] = 1
                df[i, :renlfail] = 1
            end
        end

        # Exclusions according to DRG
        if drg != nothing
            drgflg = drgmap[df[i, drg]]
            if drgflg == 1 # CARDDRG
                df[i, [:chf, :valve, :pulmcirc, :perivasc]] .= 0
                htn = 0

                # if htncx was coded with any of these tmpflg values, reset it to 0
                if any(x -> tmpflg[x], [2, 3, 6, 7, 8, 9, 10])
                    htncx = 0
                end

                # if chf was coded with tmpflg == 3
                if any(x -> tmpflg[x], [3, 7, 9])
                    df[i, :chf] = 0
                end
            elseif drgflg == 2 # PERIDRG
                df[i, :perivasc] = 0
            elseif drgflg == 3 # RENALDRG
                if any(x -> tmpflg[x], [4, 6, 7, 10])
                    htncx = 0
                end
                if any(x -> tmpflg[x], [5, 8, 9])
                    htncx = 0
                    df[i, :renlfail] = 0
                end
            elseif drgflg == 4 # NERVDRG
                df[i, :neuro] = 0
            elseif drgflg == 5 # CEREDRG
                df[i, :para] = 0
            elseif drgflag == 6 # PULMDRG
                df[i, :chrnlung] = 0
                df[i, :pulmcirc] = 0
            elseif drgflag == 7 # DIABDRG
                df[i, :dm] = 0
                df[i, :dmcx] = 0
            elseif drgflag == 8 # HYPODRG
                df[i, :hypothy] = 0
            elseif drgflag == 9 # RENFDRG
                df[i, :renlfail] = 0
            elseif drgflg == 10 # LIVERDRG
                df[i, :liver] = 0
            elseif drgflg == 11 # ULCEDRG
                df[i, :ulcer] = 0
            elseif drgflg == 12 # HIVDRG
                df[i, :aids] = 0
            elseif drgflg == 13 # LEUKDRG
                df[i, :lymph] = 0
            elseif drgflg == 14 # CANCDRG
                df[i, :tumor] = 0
            elseif drgflg == 15 # ARTHDRG
                df[i, :arth] = 0
            elseif drgflg == 16 # NUTRDRG
                df[i, :wghtloss] = 0
                df[i, :obese] = 0
                df[i, :wghtloss] = 0
                df[i, :lytes] = 0
            elseif drgflg == 17 # ANEMDRG
                df[i, :bldloss] = 0
                df[i, :anemdef] = 0
            elseif drgflg == 18 # ALCDRG
                df[i, :alcohol] = 0
                df[i, :drug] = 0
            elseif drgflg == 20 # HTNCXDRG
                htncx = 0
            elseif drgflg == 22 # PSYDRG
                df[i, :psych] = 0
            elseif drgflg == 23 # OBESEDRG
                df[i, :obese] = 0
            elseif drgflg == 24 # DEPRSDRG
                df[i, :depress] = 0
            end
        end

        # mutually exclusive conditions
        if htncx == 1
            htn = 0
        end
        if df[i, :mets] == 1
            df[i, :tumor] = 0
        end
        if df[i, :dmcx] == 1
            df[i, :dm] = 0
        end

        # htn_c
        # if htn == 1 || htncx == 1
        #     df[i, :htnc] = 1
        # end
    end
end



"""
    elixhauser10!(df, icdvars::Vector; poa = [], icdver = nothing)

Produces 39 variables of Int8 type that identifies comorbidities in the
Elixhauser Comorbidity Index. 

- `icdvars` are a vector of variable names that contain ICD-10 diagnostic codes. 

- `poa` indicates whether a condition was present on admission. If specified, these variables must
match the number of `icdvars`. This program assumes that the order of `icdvars` and `poa`
variables are identical so that the first ICD variable uses the first POA indicator, etc.

- `icdver` is the indicator variable of the ICD-10 version. This program supports version 33 through 43.
This variable should contain numeric version number (33 - 43). If not specified, the program will assume version 43.
You can create the ICD-10 version based on the date of the record using the `icd10version` function.
"""
function elixhauser10!(df, icdvars::Vector; poa = [], icdver = nothing)
    POA = length(poa) == length(icdvars) ? true : false

    # create a POA vector
    if POA
        poavars = falses(Int8, length(poa))
        cmr_cbvd_npoa = 0
        cmr_cbvd = 0
    end

    # load ICD-10 data
    elixdata = load(joinpath(@__DIR__,"..","data", "elixhauser_v10.jld2"))
    # dd = elixdata["dd"] # ICD to disease mapping
    # condnm = elixdata["conddesc"] # condition names
    # description = elixdata["desc"] # condition descriptions
    # poaexempt = elixdata["poaexempt"] # 20 POA exempt conditions (1,2,4,6,7,8,9,10,14,15,16,17,18,20,21,24,28,30,35,36)
    if icdver == nothing
        poaxmpt_codes = elixdata["poaxmpt_codes"]["v43"] # POA exempt ICD-10 codes
    end
    for v in condnm
        df[:, v] = zeros(Int8, nrow(df))
        label!(df, v, dondesc[v])
    end

    # iterate over all source records
    for i in 1:nrow(df)
        if idver != nothing && df[i,idver] != 43
            poaxmpt_codes = elixdata["poaxmpt_codes"][string("v", df[i,idver])]
        end
        if POA
            poavars = [ !ismissing(x) || in(x, ["Y","W",1,true]) ? true : false for x in df[i,poa]] 
            cmr_cbvd_npoa = 0
            cmr_cbvd = 0
        end
        for (k,icd) in enumerate(df[i,icdvars])
            if ismissing(icd) || icd in (""," ")
                continue
            end
            if haskey(dd, icd)
                # find the index for the ICD-10 code
                idx = dd[icd] # an ICD code can be mapped to 2 conditions so we need to iterate the returned tuple
                for j in idx
                    # if POA is not specified (e.g., outpatient data do not have POA codes)
                    # or if the condition is POA exempt (20 conditions are POA exempt)
                    # or if the ICD-10 code is POA exempt (255 codes are POA exempt)
                    # or if the ICD-10 code has a matching POA code showing the condition was present on admission
                    if j < 39 && any(POA, poaexempt[j], in(icd,poaxmpt_codes), poavars[k])
                        vv = condnm[j]
                        df[i, vv] = 1
                    else
                        # these are 

                    end
                end
            end
        end
        # mutually exclusive conditions (CBVD is not coded)
        if df[i, :diab_cx] == 1
            df[i, :diabimcx] = 0
        end
        if df[i, :htn_cx] == 1
            df[i, :htn_uncx] = 0
        end
        if df[i,:cancer_mets] == 1
            df[i,:cancer_solid] = 0
            df[i,:cancer_nsitu] = 0
        end
        if df[i,:cancer_solid] == 1
            df[i,:cancer_nsitu] = 0
        end
        if df[i, :liver_sev] == 1
            df[i,:liver_mld] = 0
        end
        if df[i,:renlfl_sev] == 1
            df[i,:renlfl_mod] = 0
        end
    end
end

"""
    icd10version(rdates::Vector{Dates.Date})

Returns ICD-10 version number inferred based on the servicde date.
"""
function icd10version(rdates::Vector{Dates.Date})
    rvec = Vector{Union{Missing,Int8}}(missinng, length(rdates))
    for i in 1:length(rdates)
        if ismissing(rdates[i])
            rvec[i] = missing
            continue
        end
        t = year(rdates[i]) * 10 + quarterofyear(rdates[i])
        if t in (20154, 20161, 20162, 20163)
            rvec[i] = 33
        elseif t in (20164, 20171, 20172, 20173)
            rvec[i] = 34
        elseif t in (20174, 20181, 20182, 20183)
            rvec[i] = 35
        elseif t in (20184, 20191, 20192, 20193)
            rvec[i] = 36
        elseif t in (20194, 20201, 20202, 20203)
            rvec[i] = 37
        elseif t in (20204, 20211, 20212, 20213)
            rvec[i] = 38
        elseif t in (20214, 20221, 20222, 20223)
            rvec[i] = 39
        elseif t in (20224, 20231, 20232, 20233)
            rvec[i] = 40
        elseif t in (20234, 20241, 20242, 20243)
            rvec[i] = 41
        elseif t in (20244, 20251, 20252, 20253)
            rvec[i] = 42
        elseif t in (20254, 20261, 20262, 20263)
            rvec[i] = 43
        end
    end
    return rvec
end

function elixhauser_combine(df, id, condvars::Vector)

    # create output dataset
    uniqid = sort(unique(df[:, id]))
    odf = DataFrame(id=uniqid)
    rename!(odf, id => id)

    # create comorbidity variables
    for v in condvars
        odf[i, v] = zeros(Int8, nrow(odf))
    end

    # summarize the patient records
    # DO NOT SORT. FOR LARGE DATASETS, IT MAY TAKE TOO MUCH TIME
    for subdf in groupby(df, id)

        pid = subdf[1, id]
        i = findfirst(x -> x == pid, uniqid)
        for v in condvars
            odf[i, v] = maximum(skipmissing(subdf[:, v]))
        end
    end

    # attach variable labels
    labels = Dict(names(df) .=> labels(df))
    for v in names(odf)
        if haskey(labels, v)
            label!(odf, v, labels[v])
        end
    end

    return odf
end

