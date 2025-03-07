#--------------------------------------------------------------------------
# t-test
#--------------------------------------------------------------------------
struct TTEST
    title::String
    colnms::Vector{String}
    array::Array
    levels::Vector
    μ0::Real
    t::Float64
    dof::Float64 # for UnequalVarianceTTest
    p_left::Float64
    p_twosided::Float64
    p_right::Float64
    paired::Bool
    welch::Bool
    sig::Float64
    varname::String
end

"""
    ttest(df::DataFrame,varname::Symbol;byvar::Union{Nothing,Symbol}=nothing,sig=95,welch=false)
    ttest(df::DataFrame,varname1::Symbol,varname2::Symbol;paired::Bool=false,welch=false,sig=95)
    ttest(df::DataFrame,varname::Symbol,μ0::Real;sig=95)
    ttest(var1::AbstractVector,var2::AbstractVector;paired::Bood=false,welch::Bool=false,sig=95)

Performs one-sample and two-sample t tests and produces a summary table with t statistic,
degree of freedom, and associated p-values. For paired t tests, specify `paired = true`.
For two samples with unequal variances, use `welch = true` option.

### Example 1: One-sample t test comparing a variable to a value (μ under H₀)

```
julia> t = ttest(auto, :price, 25000)
One-sample t test
| Variable |      N |      Mean |        SD |       SE |    95% LB |    95% UB |
├──────────┼────────┼───────────┼───────────┼──────────┼───────────┼───────────┤
|    price |     74 | 6165.2568 | 2949.4959 | 342.8719 | 5481.9140 | 6848.5995 |
diff = mean(price)
H₀: diff = 25000
t = -54.9323 (df = 73)

Hₐ: diff < 25000                Hₐ: diff != 25000               Hₐ: diff > 25000
Pr(T < t) = 0.0000           Pr(|T| > |t|) = 0.0000           Pr(T > t) = 1.0000
```

### Example 2: Two-sample t test by groups (equal variances)

```
julia> ttest(auto, :price, byvar = :foreign)
Two-sample t test with equal variances
|  foreign |      N |      Mean |        SD |       SE |    95% LB |    95% UB |
├──────────┼────────┼───────────┼───────────┼──────────┼───────────┼───────────┤
|        0 |     52 | 6072.4231 | 3097.1043 | 429.4911 | 5210.1837 | 6934.6624 |
|        1 |     22 | 6384.6818 | 2621.9151 | 558.9942 | 5222.1898 | 7547.1738 |
├──────────┼────────┼───────────┼───────────┼──────────┼───────────┼───────────┤
| combined |     74 | 6165.2568 | 2949.4959 | 342.8719 | 5481.9140 | 6848.5995 |
├──────────┼────────┼───────────┼───────────┼──────────┼───────────┼───────────┤
|     diff |        | -312.2587 |           | 754.4488 | -1816.225 | 1191.7075 |

diff = mean(0) - mean(1)
H₀: diff = 0
t = -0.4139 (df = 72)

Hₐ: diff < 0                      Hₐ: diff != 0                     Hₐ: diff > 0
Pr(T < t) = 0.3401           Pr(|T| > |t|) = 0.6802           Pr(T > t) = 0.6599
```

### Example 3: Two-sample t test with unequal variances

```
julia> ttest(auto,:price,by=:foreign, welch=true)
Two-sample t test with unequal variances
|  foreign |      N |      Mean |        SD |       SE |    95% LB |    95% UB |
├──────────┼────────┼───────────┼───────────┼──────────┼───────────┼───────────┤
|        0 |     52 | 6072.4231 | 3097.1043 | 429.4911 | 5210.1837 | 6934.6624 |
|        1 |     22 | 6384.6818 | 2621.9151 | 558.9942 | 5222.1898 | 7547.1738 |
├──────────┼────────┼───────────┼───────────┼──────────┼───────────┼───────────┤
| combined |     74 | 6165.2568 | 2949.4959 | 342.8719 | 5481.9140 | 6848.5995 |
├──────────┼────────┼───────────┼───────────┼──────────┼───────────┼───────────┤
|     diff |        | -312.2587 |           | 704.9376 | -1730.856 | 1106.3386 |

diff = mean(0) - mean(1)
H₀: diff = 0
t = -0.4430 (df = 46.4471)

Hₐ: diff < 0                      Hₐ: diff != 0                     Hₐ: diff > 0
Pr(T < t) = 0.3299           Pr(|T| > |t|) = 0.6599           Pr(T > t) = 0.6701
```

### Example 4: Two-sample paired t test comparing two variables

```
julia> ttest(fuel, :mpg1, :mpg2, paired = true)
Paired t test
| Variable |      N |     Mean |       SD |       SE |   95% LB |   95% UB |
├──────────┼────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
|     mpg1 |     12 | 21.00000 |  2.73030 |  0.78817 | 19.26525 | 22.73475 |
|     mpg2 |     12 | 22.75000 |  3.25087 |  0.93845 | 20.68449 | 24.81551 |
├──────────┼────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
|     diff |     12 | -1.75000 |  2.70101 |  0.77971 | -3.46614 | -0.03386 |

diff = mean(mpg1) - mean(mpg2)
H₀: diff = 0
t = -2.2444 (df = 11)

Hₐ: diff < 0                   Hₐ: diff != 0                    Hₐ: diff > 0
Pr(T < t) = 0.0232         Pr(|T| > |t|) = 0.0463         Pr(T > t) = 0.9768
```

### Example 5: Two-sample unpaired t test comparing two variables

```
julia> ttest(fuel, :mpg1, :mpg2)
Two-sample t test with equal variances
| Variable |      N |     Mean |       SD |       SE |   95% LB |   95% UB |
├──────────┼────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
|     mpg1 |     12 | 21.00000 |  2.73030 |  0.78817 | 19.26525 | 22.73475 |
|     mpg2 |     12 | 22.75000 |  3.25087 |  0.93845 | 20.68449 | 24.81551 |
├──────────┼────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
| combined |     24 | 21.87500 |  3.06895 |  0.62645 | 20.57909 | 23.17091 |
├──────────┼────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
|     diff |        | -1.75000 |          |  1.22552 | -4.29157 |  0.79157 |

diff = mean(mpg1) - mean(mpg2)
H₀: diff = 0
t = -1.4280 (df = 22)

Hₐ: diff < 0                   Hₐ: diff != 0                    Hₐ: diff > 0
Pr(T < t) = 0.0837         Pr(|T| > |t|) = 0.1673         Pr(T > t) = 0.9163
```

"""
function ttest(df::DataFrame, varname::Symbol; by::Symbol = nothing, sig = 95, welch = false)
    if by == nothing
        error("`by` is required.")
    end

    dft = dropmissing( df[!, [varname, by]] )

    lev = sort(unique(dft[!,by]))
    if length(lev) != 2
        error(by," must have two levels; it has ",length(lev), " levels.")
    end

    val1 = dft[dft[!, by] .== lev[1],varname]
    val2 = dft[dft[!, by] .== lev[2],varname]
    if length(val1) == 0 || length(val2) == 0
        error("On or both groups have zero observations")
    end
    return ttest(val1,val2,paired=false,welch=welch,sig=sig,levels=lev, by = by)
end
function ttest(df::DataFrame, var1::Symbol, var2::Symbol; sig = 95, paired = false, welch = false)

    if paired == true
        df2 = DataFrames.dropmissing(df[:,[var1,var2]])
        x = df2[:, var1]
        y = df2[:, var2]
    else
        x = collect(skipmissing(df[!,var1]))
        y = collect(skipmissing(df[!,var2]))
    end

    length(x) > 0 && length(y) > 0 || error("One or both variables are empty")

    return ttest(x, y, paired=paired, welch=welch, levels=[var1,var2], sig = sig, by = nothing)
end
function ttest(x::AbstractVector, y::AbstractVector; 
    paired::Bool=false,
    welch::Bool=false,
    sig=95,
    levels=Any[:x,:y],
    by = by)

    paired && length(x) != length(y) && error("Paired t test requires equal lengths in input vectors")

    if paired == true
        title = "Paired t test"
        tt = OneSampleTTest(x, y)
    elseif welch == false
        title = "Two-sample t test with equal variances"
        tt = EqualVarianceTTest(x, y)
    else
        title = "Two-sample t test with unequal variances"
        tt = UnequalVarianceTTest(x, y)
    end

    # compute standard errors and confidence intervals
    N = Vector{Union{Missing,Int64}}(missing, 4)
    MEAN = Vector{Float64}(undef, 4)
    SD = Vector{Union{Missing,Float64}}(missing, 4)
    SE = Vector{Float64}(undef, 4)
    LB = Vector{Float64}(undef, 4)
    UB = Vector{Float64}(undef, 4)

    val = [x,y]
    for i = 1:2
        N[i] = length(val[i])
        MEAN[i] = mean(val[i])
        SD[i] = std(val[i])
        SE[i] = SD[i] / sqrt(N[i])
        LB[i],UB[i] = StatsAPI.confint(OneSampleTTest(val[i]))
    end

    # combined
    N[3] = N[1] + N[2]
    MEAN[3] = mean(vcat(val[1],val[2]))
    SD[3] = std(vcat(val[1],val[2]))
    SE[3] = SD[3] / sqrt(N[3])
    LB[3],UB[3] = StatsAPI.confint(OneSampleTTest(vcat(val[1],val[2])))

    # diff
    N[4] = paired ? N[1] : missing
    MEAN[4] = tt.xbar
    SD[4] = paired ? tt.stderr*sqrt(N[4]) : missing
    SE[4] = tt.stderr
    LB[4],UB[4] = StatsAPI.confint(tt)

    return TTEST(title,
            ["N", "Mean", "SD", "SE", string(sig,"% LB"), string(sig,"% UB")],
            [N MEAN SD SE LB UB],
            levels,
            tt.μ0,
            tt.t,
            tt.df,
            pvalue(tt, tail = :left),
            pvalue(tt),
            pvalue(tt, tail = :right),
            paired,
            welch,
            sig,
            string(by))
end
function ttest(df::DataFrame, varname::Symbol, μ0::Real; sig=95)
    v = collect(skipmissing(df[!,varname]))
    ttest(v, μ0, varname = varname, sig=sig)
end
function ttest(var::AbstractVector, μ0::Real = 0; varname = nothing, sig = 95)

    N = length(var)
    if N == 0
        error("Vector is empty")
    end
    MEAN = mean(var)
    SD = std(var)

    tt = OneSampleTTest(var, μ0)
    title = "One-sample t test"

    # compute standard errors and confidence intervals
    SE = SD / sqrt(N)
    (LB, UB) = StatsAPI.confint(tt)

    return TTEST(title,
        ["N", "Mean", "SD", "SE", string(sig,"% LB"), string(sig,"% UB")],
        [[N] [MEAN] [SD] [SE] [LB] [UB]],
        [],
        tt.μ0,
        tt.t,
        tt.df,
        pvalue(tt, tail = :left),
        pvalue(tt),
        pvalue(tt, tail = :right),
        false,
        false,
        sig,
        varname == nothing ? "x" : string(varname))

end
function Base.show(io::IO, t::TTEST)
    println(io,"\n",t.title,"\n")

    if t.title == "One-sample t test"
        pretty_table(io, 
            t.array,
            header=t.colnms,
            row_labels= [t.varname == nothing ? "" : string(t.varname)],
            row_label_column_title="Variable",
            formatters=(ft_printf("%.5f", [2, 3, 4, 5, 6]), ft_printf("%.0f", [1])),
            hlines=[0, 1, 2],
            vlines=[1]
        )
        println(io, "\nmean = mean(", t.varname, ")")
        println(io, "H₀: mean = ", t.μ0)
        println(io, "t = ", @sprintf("%.5f", t.t), " (df = ", @sprintf("%.0f", t.dof), ")\n")
    else
        pretty_table(io,
            t.array,
            header=t.colnms,
            row_labels = vcat(t.levels,"combined","diff"),
            row_label_column_title = t.varname == nothing ? "Variable" : string(t.varname),
            formatters = (ft_printf("%.4f",[2,3,4,5,6]),ft_printf("%.0f",[1]), ft_nomissing),
            hlines = [1,3,4],
            vlines = [1]
        )
        println(io, "\ndiff = mean(", t.levels[1], ") - mean(", t.levels[2], ")")
        println(io, "H₀: diff = 0")
        println(io, "t = ", @sprintf("%.5f", t.t), " (df = ", @sprintf("%.5f",t.dof), ")\n")
    end

    pretty_table(io,
        [string("P = ", @sprintf("%.4f",t.p_left)) string("P = ", @sprintf("%.4f",t.p_twosided)) string("P = ", @sprintf("%.4f",t.p_right))],
        header = ["Hₐ: diff < 0       ","       Hₐ: diff != 0       ","       Hₐ: diff > 0"],
        # formatters = (ft_printf("%.4f")),
        alignment = [:l,:c,:r],
        hlines = :none,
        vlines = :none
    )

end
