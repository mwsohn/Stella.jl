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
    kaplanmeier(df, event; by = nothing)

Plots Kaplan-Meier estimates.
"""
function kaplanmeier(df, event, by=nothing, kwargs...)

    plt = nothing
    if by == nothing
        km = fit(KaplanMeier, df[:, event])
        plot(vcat(0, km.events.time), vcat(1.0, km.survival), linetype=:steppost, ylims=(0, 1), kwargs...)
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
                    ylabel="Survival estimate",
                    label=string(v),
                    kwargs...)

            else
                Plots.plot!(plt, vcat(0, kvec[i].events.time),
                    vcat(1.0, kvec[i].survival),
                    linetype=:steppost,
                    ylims=(0, 1),
                    label=string(v),
                    kwargs...)

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