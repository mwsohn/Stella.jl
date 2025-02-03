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
    if interaction
        fm = @eval @formula($dep ~ 1 + $cat1 + $cat2 + $cat1 * $cat2)
    else
        fm = @eval @formula($dep ~ 1 + $cat1 + $cat2)
    end
    return anova(_df, fm, type = type)
end
function anova(_df::AbstractDataFrame, fm; type = 1)
    dep = fm.lhs.sym
    intercept = isa(fm.rhs[1], ConstantTerm) && fm.rhs[1].n == 1 ? true : false
    cats = Vector{Symbol}()
    nlev = Vector{Int}()
    interaction = false
    for i = 1:length(fm.rhs)
        if i == 1 && intercept
            continue
        end
        if length(terms(fm.rhs[i])) == 1
            push!(cats,fm.rhs[i].sym)
            isa(_df[:, cats[1]], CategoricalArray) || throw(ArgumentError(cats[i], " must be a Categorical Array"))
            push!(nlev,length(unique(skipmissing(_df[:,cats[i]]))))
        elseif isa(fm.rhs[i], InteractionTerm) && length(terms(fm.rhs[i])) == 2
            interaction = true
            push!(nlev,(nlev[1]-1)*(nlev[2]-1)+1)
        end
    end
    ba = completecases(_df[:,vcat(dep, cats)])
    df2 = _df[ba,vcat(dep, cats)]
    mm = modelmatrix(fm,df2)
    X = hcat(mm, df2[:,dep])
    if intercept == false
        X = hcat(ones(Float64,size(X,1)),X)
    end
    XX = X'X
    SS = type == 1 ? SSTypeI(XX, nlev) : SSTypeII(XX, nlev)
    tdf = nrow(df2) - 1
    mdf = sum(nlev) - length(nlev)
    DF = vcat(mdf, nlev .- 1, tdf - mdf, tdf )
    rdf = tdf - mdf
    MSS = SS ./ DF
    rms = MSS[end-1]
    
    Source = interaction ? 
        ["Model", string(cats[1]), string(cats[2]), string(cats[1], " & ", cats[2]), "Residual", "Total"] :
        ["Model", string(cats[1]), string(cats[2]), "Residual", "Total"]

    return AOV(
        type == 1 ? "Type I" : "Type II",
        Source,
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
function SSTypeII(XX,nlev)
    A = copy(XX)
    (r,c) = size(A)
    n = length(nlev)
    SS = zeros(Float64,n+3)
    sweep!(A,1)
    # TSS
    SS[n+3] = copy(A[r,c])
    # RSS
    sweep!(A,2:n)
    SS[n+2] = copy(A[r,c])
    # invert factors
    pos = 2
    for (i,v) in enumerate(nlev)
        B = copy(A)
        sweep!(B,pos:(pos+v-2),true)
        pos += (v-1)
        SS[i+1] = B[r,c] - SS[n+2]
    end
    sweep!(A,2:n,true)
    # MSS
    SS[1] = A[r,c] - SS[n+2]
    return SS
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
