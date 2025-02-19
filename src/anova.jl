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

## Examples
We will sue `auto` dataset downloaded from https://vincentarelbundock.github.io/Rdatasets/csv/causaldata/auto.csv is used in this example. We want to
examine whether car prices are different by domestic- vs foreign-made or by MPG tertiles.

We will create MPG tertiles based on auto.mpg and convert them into CategoricalArrays. `xtile` and `values!` functions are both part of the "unpublished" Stella.jl package.

```
julia> auto.mpg3 = xtile(auto, :mpg, nq = 3);

julia> values!(auto,:foreign, Dict(0 => "Domestic", 1 => "Foreign"));

julia> tab(auto, :foreign)
──────────┬───────────────────────────
  foreign │ Counts   Percent  Cum Pct 
──────────┼───────────────────────────
 Domestic │     52   70.2703  70.2703
  Foreign │     22   29.7297    100.0
──────────┼───────────────────────────
    Total │     74     100.0    100.0
──────────┴───────────────────────────

julia> values!(auto,:mpg3, Dict(1 => "MPG Tertile 1", 2 => "MPG Tertile 2", 3 => "MPG Tertile 3"));

julia> tab(auto, :mpg3)
───────────────┬───────────────────────────
          mpg3 │ Counts   Percent  Cum Pct 
───────────────┼───────────────────────────
 MPG Tertile 1 │     27   36.4865  36.4865
 MPG Tertile 2 │     24   32.4324  68.9189
 MPG Tertile 3 │     23   31.0811    100.0
───────────────┼───────────────────────────
         Total │     74     100.0    100.0
───────────────┴───────────────────────────
```

### 1. Oneway ANOVA

```
julia> aov = anova(auto, :price, :mpg3)

Analysis of Variance (One-Way)

   Source │            SS  DF            MS      F       P 
──────────┼────────────────────────────────────────────────
    Model │ 136763694.954   2  68381847.477  9.743  0.0002
     mpg3 │ 136763694.954   2  68381847.477  9.743  0.0002
 Residual │ 498301701.168  71   7018333.819
──────────┼────────────────────────────────────────────────
    Total │ 635065396.122  73   8699525.974

julia> aov.pvalue[1]
0.00018235773127089433
```
### 2. Twoway ANOVA

#### 2.1. Type I Sums of Squares, No interaction
```
julia> aov = anova(auto, :price, :foreign, :mpg3, type = 1)

Analysis of Variance (Type I)

   Source │            SS  DF            MS       F         P 
──────────┼───────────────────────────────────────────────────
    Model │ 154897061.512   3  51632353.837   7.527    0.0002
  foreign │   1507382.657   1   1507382.657   0.220    0.6407
     mpg3 │ 153389678.855   2  76694839.428  11.181  < 0.0001
 Residual │ 480168334.610  70   6859547.637
──────────┼───────────────────────────────────────────────────
    Total │ 635065396.122  73   8699525.974

julia> aov.pvalue[3] # P-value for mpg3
6.112930853953449e-5
```

#### 2.2. Type II Sums of Squares, No interaction
```
julia> aov = anova(auto, :price, :foreign, :mpg3, type = 2)

Analysis of Variance (Type II)

   Source │            SS  DF            MS       F         P 
──────────┼───────────────────────────────────────────────────
    Model │ 154897061.512   3  51632353.837   7.527    0.0002
  foreign │  18133366.558   1  18133366.558   2.644    0.1085
     mpg3 │ 153389678.855   2  76694839.428  11.181  < 0.0001
 Residual │ 480168334.610  70   6859547.637
──────────┼───────────────────────────────────────────────────
    Total │ 635065396.122  73   8699525.974

julia> aov.pvalue[3] # P-value for mpg3
6.112930853953449e-5
```

#### 2.3. Type I Sums of Squares, With Interaction
```
julia> aov = anova(auto, :price, :foreign, :mpg3, type = 1, interaction = true)

Analysis of Variance (Type I)

         Source │            SS  DF            MS       F         P 
────────────────┼───────────────────────────────────────────────────
          Model │ 156625170.566   5  31325034.113   4.452    0.0014
        foreign │   1507382.657   1   1507382.657   0.214    0.6449
           mpg3 │ 153389678.855   2  76694839.428  10.901  < 0.0001
 foreign & mpg3 │   1728109.054   2    864054.527   0.123    0.8846
       Residual │ 478440225.555  68   7035885.670
────────────────┼───────────────────────────────────────────────────
          Total │ 635065396.122  73   8699525.974

julia> aov.pvalue[3] # P-value for mpg3
7.829523220630382e-5
```

#### 2.4. Type II Sums of Squares, With Interaction
```
julia> aov = anova(auto, :price, :foreign, :mpg3, type = 2, interaction = true)

Analysis of Variance (Type II)

         Source │            SS  DF            MS       F         P 
────────────────┼───────────────────────────────────────────────────
          Model │ 156625170.566   5  31325034.113   4.452    0.0014
        foreign │  18133366.558   1  18133366.558   2.577    0.1130
           mpg3 │ 153389678.855   2  76694839.428  10.901  < 0.0001
 foreign & mpg3 │   1728109.054   2    864054.527   0.123    0.8846
       Residual │ 478440225.555  68   7035885.670
────────────────┼───────────────────────────────────────────────────
          Total │ 635065396.122  73   8699525.974

julia> aov.pvalue[3] # P-value for mpg3
7.829523220630382e-5
```

#### 2.5. Type III Sums of Squares, With Interaction
```
julia> aov = anova(auto, :price, :foreign, :mpg3, type = 3, interaction = true)

Analysis of Variance (Type III)

         Source │            SS  DF            MS      F       P 
────────────────┼────────────────────────────────────────────────
          Model │ 156625170.566   5  31325034.113  4.452  0.0014
        foreign │  19148362.309   1  19148362.309  2.722  0.1036
           mpg3 │ 129967822.960   2  64983911.480  9.236  0.0003
 foreign & mpg3 │   1728109.054   2    864054.527  0.123  0.8846
       Residual │ 478440225.555  68   7035885.670
────────────────┼────────────────────────────────────────────────
          Total │ 635065396.122  73   8699525.974

julia> aov.pvalue[3] # P-value for mpg3
0.0002828218888318445
```

#### 2.6. Type III Sums of Squares, With Interaction using Formula
```
julia> aov = anova(auto, @formula(price ~ foreign + mpg3 + foreign & mpg3), type=3)

Analysis of Variance (Type III)

         Source │            SS  DF            MS      F       P 
────────────────┼────────────────────────────────────────────────
          Model │ 156625170.566   5  31325034.113  4.452  0.0014
        foreign │  19148362.309   1  19148362.309  2.722  0.1036
           mpg3 │ 129967822.960   2  64983911.480  9.236  0.0003
 foreign & mpg3 │   1728109.054   2    864054.527  0.123  0.8846
       Residual │ 478440225.555  68   7035885.670
────────────────┼────────────────────────────────────────────────
          Total │ 635065396.122  73   8699525.974

julia> aov.pvalue[3] # P-value for mpg3
0.0002828218888318445
```

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
    sweep!(A,2:(len[1]-1))
    RSS = copy(A[len...])
    MSS = TSS - RSS
    mdf = len[1] - 2
    tdf = nrow(df2) - 1
    mms = MSS/mdf
    rms = RSS / (tdf - mdf)
    tms = TSS / tdf
    pval = ccdf(FDist(mdf, tdf - mdf), mms / rms)

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
    MF = ModelFrame(fm, _df, contrasts=Dict(:foreign => EffectsCoding(), :mpg3 => EffectsCoding()))
    terms = MF.f.rhs.terms
    cats = Vector{Symbol}()
    nlev = Vector{Int}()
    for i = 2:length(terms)
        if isdefined(terms[i], :sym) && isdefined(terms[i], :contrasts)
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
