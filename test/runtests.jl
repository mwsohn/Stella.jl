

using DataFrames, Stella, RDatasets, CategoricalArrays

# ------------------------------------------------------
# t-test
# ------------------------------------------------------

@testset "TTEST" begin

    # read_stata
    auto = read_stata(download("http://www.stata-press.com/data/r13/auto.dta"))

    # one-sample ttest comparing a variable with a value
    t = ttest(auto, :price, 25000)
    @test t.t ≈ -54.9323 atol=0.001
    @test t.p_twosided ≈ 0.0000 atol=0.001

    # two-sample ttest by groups (equal variance)
    t = ttest(auto, :price, by=:foreign)
    @test t.t = -0.4139 atol=0.001
    @test t.p_twosided ≈ 0.6802 atol=0.001

    # Two-sample t test with unequal variances
    t = ttest(auto, :price, by=:foreign, welch=true)
    @test t.t ≈ -0.4430 atol=0.001
    @test t.p_twosided ≈ 0.6598 atol=0.001 

    # Two-sample paired t test comparing two variables
    fuel = read_stata(download("https://www.stata-press.com/data/r18/fuel.dta"))
    t = ttest(fuel, :mpg1, :mpg2, paired=true)
    @test t.t ≈ -2.2444 atol=0.001
    @test t.p_twosided ≈ 0.0463 atol=0.001

    # Two-sample unpaired t test comparing two variables - paired = false is the default
    t = ttest(fuel, :mpg1, :mpg2, paired=false)
    @test t.t ≈ -1.4280 atol=0.001
    @test t.p_twosided ≈ 0.1673 atol=0.001

end

# ------------------------------------------------------
# ANOVA
# ------------------------------------------------------

    