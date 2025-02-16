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
	anova(::DataFrame, depvar::Symbol, groupvar::Symbol)
	anova(::DataFrame, depvar::Symbol, groupvar1::Symbol, groupvar2::Symbol; type = 1, interaction = false)
    anova(::DataFrame, fm::FormulaTerm; type=1)
    anova(::StatsModels.TableRegressionModel)

Performs an oneway and twoway ANOVA analysis. `depvar' is a continous variable and `groupvar1' and `groupvar2` are
group variables. Currently, `anova` supports 2 different syntax. The first syntax requires specifying variable names, ANOVA type,
and interaction directly as arguments. The second requires specifying a linear regression `formula` 
(following StatsModels.jl @formula language) and ANOVA type. 

## Options:

* `df` - input data in AbstractDataFrame
* `depvar` - dependent variable (Continous)
* `groupvar`, `groupvar1`, `groupvar2` - independent variables. They must all be CategoricalArrays.
* `fm` - formula
* `type` - Type of ANOVA sums of squares. Currently, I, II, and III are supported. Specify 1 for I, 2 for II, etc. 
    These types match the type produced in SAS GLM. Type I is the same as Stata `sequential` sums of squares and Type III
    as `partial` sums of squares.
* `interaction` - indicates whether an interaction terms should be included. In the second syntax where you specify a formula,
    an interaction is included with `groupvar1 * groupvar2` as part of the formula or with `groupvar1 & groupvar2`.

## Output:
Output is a struct whose elements are:

* type - Type of sums of squares
* title - Row titles in the ANOVA table
* ss - Sums of squares
* df - Degress of freedom
* ms - Mean sums of squares
* F - F-statistic
* pvalue - P-values
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
    isa(_df[:, cat1], CategoricalArray) || throw(ArgumentError("`cat1` must be a Categorical Array"))
    isa(_df[:, cat2], CategoricalArray) || throw(ArgumentError("`cat2` must be a Categorical Array"))
    return anova(_df, interaction ? @eval(@formula($dep ~ $cat1 + $cat2 + $cat1 * $cat2)) : @eval(@formula($dep ~ 1 + $cat1 + $cat2)), type=type)
end
function anova(_df::AbstractDataFrame, fm; type = 1)
    if type == 3
        MF = ModelFrame(fm, _df, contrasts=Dict(:foreign => EffectsCoding(), :mpg3 => EffectsCoding()))
    else
        MF = ModelFrame(fm, _df, contrasts=Dict(:foreign => EffectsCoding(), :mpg3 => EffectsCoding()))
    end

    terms = MF.f.rhs.terms
    cats = Vector{Symbol}()
    nlev = Vector{Int}()
    for i = 2:length(terms)
        if isdefined(terms[i], :sym)
            isa(_df[:, terms[i].sym], CategoricalArray) || throw(ArgumentError("requires a Categorical Array"))
            push!(cats, terms[i].sym)
            push!(nlev, length(terms[i].contrasts.levels)-1)
        else
            (sy, len) = get_sym_lev(terms[i].terms)
            push!(nlev, len)
            push!(cats, Symbol(sy))
        end
    end

    MM = modelmatrix(MF)
    X = hcat(MM, MF.data[1])
    XX = X'X
    if type == 1
        Type = "Type I"
        SS = SSTypeI(XX, nlev)
    elseif type == 2
        Type = "Type II"
        SS = SSTypeII(XX, nlev)
    elseif type == 3
        Type = "Type III"
        SS = xSSTypeIII(XX, nlev)
    end
    tdf = size(MM,1) - 1
    mdf = sum(nlev)
    DF = vcat(mdf, nlev, tdf - mdf, tdf)
    rdf = tdf - mdf
    MSS = SS ./ DF
    rms = MSS[end-1]

    return Stella.AOV(
        Type,
        vcat("Model", string.(cats), "Residual", "Total"),
        SS,
        DF,
        MSS,
        [i <= length(MSS) - 2 ? MSS[i] / rms : missing for i in 1:length(MSS)],
        [i <= length(MSS) - 2 ? ccdf(FDist(DF[i], rdf), MSS[i] / rms) : missing for i in 1:length(MSS)]
    )

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
function SSTypeI(XX, nlev)
    (r, c) = size(XX)
    n = length(nlev)
    SS = zeros(Float64, n + 3)
    # TSS
    sweep!(XX, 1)
    SS[n+3] = copy(XX[r, c])
    pos = 2
    for (i, v) in enumerate(nlev)
        sweep!(XX, pos:(pos+v-1))
        pos += v
        SS[i+1] = SS[n+3] - XX[r, c] - sum(SS[1:i])
    end
    # MSS
    SS[1] = sum(SS[2:n+1])
    # RSS
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
    sweep!(XX, 2:sum(nlev[1:2])+1)
    SS[n+1] = XX[r, c] - SS[n+2]

    # SSA and SSB
    sweep!(A, 1:sum(nlev[1:2])+1)
    rss2 = copy(A[r, c])
    pos = 2
    for (i, v) in enumerate(nlev[1:2])
        B = copy(A)
        sweep!(B, pos:(pos+v-1), true)
        pos += v
        SS[i+1] = B[r, c] - rss2
    end
    return SS
end
function xSSTypeIII(XX, nlev)
    (r, c) = size(XX)
    n = length(nlev)
    SS = zeros(Float64, n + 3)
    # TSS
    sweep!(XX, 1)
    SS[n+3] = copy(XX[r, c])
    # RSS
    sweep!(XX, 2:c-1)
    SS[n+2] = copy(XX[r, c])
    A = copy(XX)
    # MSS
    sweep!(XX, 2:c-1, true)
    SS[1] = XX[r, c] - SS[n+2]

    pos = 2
    for (i, v) in enumerate(nlev)
        B = copy(A)
        sweep!(B, pos:(pos+v-1), true)
        pos += v
        SS[i+1] = B[r, c] - SS[n+2]
    end
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
