struct AOV
    type::String
    title::Vector{String}
    ss::Vector{Float64}
    df::Vector{Int64}
    ms::Vector{Float64}
    F::Vector{Union{Missing,Float64}}
    pvalue::Vector{Union{Missing,Float64}}
end

"""
	anova(::DataFrame, contvar::Symbol, groupvar::Symbol)
    anova(::StatsModels.TableRegressionModel)

Produces an one-way ANOVA table. `contvar' is a continous variable and `groupvar' is
the group variable. 
"""
function anova(_df::AbstractDataFrame, dep::Symbol, cat::Symbol)
    isa(_df[:,cat], CategoricalArray) || throw(ArgumentError("`cat` must be a Categorical Array"))
    ba = completecases(_df[:,[dep,cat]])
    df2 = _df[ba,[dep, cat]]
    fm = @eval @formula($dep ~ 1 + $cat)
    mm = modelmatrix(fm, df2)
    X = hcat(mm,df2[!, dep])
    XX = X'X
    A = copy(XX)
    len = size(XX)
    sweep!(A,1)
    TSS = copy(A[len...])
    sweep!(A,2)
    RSS = copy(A[len...])
    MSS = TSS - RSS
    mdf = len[1] - 2
    tdf = nrow(df2) - 1
    mms = MSS/mdf
    rms = RSS / (tdf - mdf)
    tms = TSS / tdf
    pval = ccdf(FDist(tdf - mdf, mdf), mms / rms)

    return AOV(
        "One-Way",
        ["Model", string(cat), "Residual", "Total"],
        [MSS, MSS, RSS, TSS],
        [mdf, mdf, tdf - mdf, tdf],
        [mms, mms, rms, tms],
        [mms / rms, mms / rms, missing, missing ],
        [ pval, pval, missing, missing]
    );
end
function anova(_df::AbstractDataFrame, dep::Symbol, cat1::Symbol, cat2::Symbol; type = 1, interaction = false)
    return anova(_df, interaction ? @eval(@formula($dep ~ $cat1 + $cat2 + $cat1 * $cat2)) : @eval(@formula($dep ~ 1 + $cat1 + $cat2)), type=type)
end
function anova(_df::AbstractDataFrame, fm; type = 1)

    if type == 3
        MF = ModelFrame(fm, _df, contrasts = Dict(:foreign => EffectsCoding(), :mpg3 => EffectsCoding()))
    else
        MF = ModelFrame(fm, _df, contrasts=Dict(:foreign => EffectsCoding(), :mpg3 => EffectsCoding()))
    end

    rhs = MF.f.rhs.terms
    cats = Vector{Symbol}()
    nlev = Vector{Int}()
    for i = 2:length(rhs)
        if isdefined(rhs,:sym)
            push!(cats,rhs[i].sym)
            push!(nlev, length(rhs[i].contrasts.levels))
        else
            (sy, len) = get_sym_lev(rhs[i].terms)
            push!(nlev,len)
            push!(cats,sy)
        end
    end

    MM = modelmatrix(MF)
    X = hcat(MM, MM.data[1])
    XX = X'X
    if type == 1
        Type = "Type I"
        SS = SSTypeI(XX, nlev)
    elseif type == 2
        Type = "Type II"
        SS = SSTypeII(XX, nlev)
    elseif type == 3
        Type = "Type III"
        SS = SSTypeIII(XX, nlev)
    end
    tdf = nrow(df2) - 1
    mdf = sum(nlev) - length(nlev)
    DF = vcat(mdf, nlev .- 1, tdf - mdf, tdf )
    rdf = tdf - mdf
    MSS = SS ./ DF
    rms = MSS[end-1]

    return AOV(
        Type,
        vcat("Model", cats, "Residual", "Total"),
        SS,
        DF,
        MSS,
        [ i <= length(MSS) - 2 ? MSS[i] / rms : missing for i in 1:length(MSS)],
        [ i <= length(MSS) - 2 ? ccdf(FDist(DF[i], rdf), MSS[i] / rms) : missing for i in 1:length(MSS) ]
    )

end
function SSTypeI(XX,nlev)
    (r,c) = size(XX)
    n = length(nlev)
    SS = zeros(Float64,n+3)
    sweep!(XX,1)
    SS[n+3] = copy(XX[r,c])
    pos = 2
    for (i,v) in enumerate(nlev)
        sweep!(XX,pos:(pos+v-2))
        pos += (v-1)
        SS[i+1] = SS[n+3] - XX[r,c] - sum(SS[1:i])
    end
    SS[1] = sum(SS[2:n+1])
    SS[n+2] = SS[n+3] - SS[1]
    return SS
end
function SSTypeII(XX, nlev)
    (r, c) = size(XX)
    n = length(nlev)
    SS = zeros(Float64, n + 3)
    A = copy(XX)
    # TSS
    sweep!(XX, 1)
    SS[n+3] = copy(XX[r, c])
    # RSS
    sweep!(XX, 2:c-1)
    SS[n+2] = copy(XX[r, c])
    # MSS
    sweep!(XX, 2:c-1, true)
    SS[1] = XX[r, c] - SS[n+2]
    # SSAB
    sweep!(XX, 2:sum(nlev[1:2])-1)
    SS[n+1] = XX[r, c] - SS[n+2]

    # sweep again
    sweep!(A, 1:sum(nlev[1:2])-1)
    rss2 = copy(A[r, c])
    pos = 2
    for (i, v) in enumerate(nlev[1:2])
        B = copy(A)
        sweep!(B, pos:(pos+v-2), true)
        pos += (v - 1)
        SS[i+1] = B[r, c] - rss2
    end
    return SS
end
function SSTypeIII(XX, nlev)
    (r, c) = size(XX)
    n = length(nlev)
    SS = zeros(Float64, n + 3)
    A = copy(XX)
    # TSS
    sweep!(XX, 1)
    SS[n+3] = copy(XX[r, c])
    # RSS
    sweep!(XX, 2:c-1)
    SS[n+2] = copy(XX[r, c])
    # MSS
    sweep!(XX, 2:c-1, true)
    SS[1] = XX[r, c] - SS[n+2]
    # SSAB
    sweep!(XX, 2:sum(nlev[1:2])-1)
    SS[n+1] = XX[r, c] - SS[n+2]

    sweep!(A, 1:sum(nlev[1:2])-1)
    rss2 = copy(A[r, c])
    pos = 2
    for (i, v) in enumerate(nlev[1:2])
        # need to fix this part
        B = copy(A)
        sweep!(B, pos:(pos+v-2), true)
        pos += (v - 1)
        SS[i+1] = B[r, c] - rss2
    end
    return SS
end
function get_sym_lev(terms)
    vstr = ""
    levs = 1
    for i = 1:length(terms)
        if isdefined(terms[i], :sym)
            vstr = string(vstr, i > 1 ? " & " : "", terms[i].sym)
            levs = levs * (length(terms[i].contrasts.levels) - 1)
        end
    end
    return vstr, levs
end

function anova(glmmodel)
    tss = nulldeviance(glmmodel)
    rss = deviance(glmmodel)
    mss = tss - rss
    tdf = nobs(glmmodel) - 1
    mdf = GLM.dof(glmmodel) - 1
    rdf = tdf - mdf
    F = (mss / mdf) / (rss / rdf)
    return AOV(
        "Regression Model",
        ["Model", "Residual", "Total"],
        [mss, rss, tss],
        [mdf, rdf, tdf],
        [(tss - rss) / mdf, rss / rdf, tss / tdf],
        [F, missing, missing],
        [Distributions.ccdf(Distributions.FDist(mdf,rdf), F), missing, missing]
    );
end
function Base.show(io::IO, a::AOV)
    n = length(a.ss)
    pstr = [ x < 0.0001 ? "< 0.0001" : @sprintf("%.4f", x) for x in skipmissing(a.pvalue) ]
    println(io, "\nAnalysis of Variance (",a.type,")\n")
    pretty_table(io, 
            DataFrame(
            Source=a.title, 
            SS=a.ss, 
            DF=a.df, 
            MS=a.ms,
            F=a.F,
            P=vcat(pstr, missing, missing)
        );
        formatters=(ft_nomissing, ft_printf("%.3f", [2, 4, 5])),
        hlines=[1, n],
        vlines=[1],
        show_subheader=false
    )
end
