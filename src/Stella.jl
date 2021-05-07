VERSION >= v"1.0.0"

__precompile__()

module Stella

################################################################################
##
## Dependencies
##
################################################################################

using DataFrames, Distributions, StatsBase, GLM, NamedArrays, HypothesisTests,
        DataStructures, FreqTables, ReadStat, Labels, LinearAlgebra, Printf, 
        Glob, Dates, Arrow, CategoricalArrays

##############################################################################
##
## Exported methods and types (in addition to everything reexported above)
##
##############################################################################

export  read_stata,  # read stata 13 and 14 files into DF
        read_stata!, # read stata 13 and 14 files into DF
        dfcompress!, # compress DF
        dacompress,  # compress a DataArray
        desc,        # list variables with Stata labels and value labels
        nmissing,    # returns the number of missing values in an AbstractArray
        dflist,      # list
        dfmerge,     # merge two dataframes
        dfsample,    # select a sample from a df
        univariate,  # univariate statistics in DF
        tabstat,     # unviariate statistics by subgroups
        substat,     # attach univariate stat to original DF by subgroups
        xtile,       # create variable that classify a column into percentiles
        rowsum,      # create variable that sums rows
        strval,      # convert floats or ints to strings
        chi2test,    # compute chisquare statistics from na.array
        prepend_spaces, append_spaces, # create fixed length strings
        smallest, largest, # list smallest and largest values in a DA
        pickone,     # create a binary variable that identifies one records in a subgroup
        p5, p10, p25, p50, p75, p90, p95, # percentiles
        cv,          # coefficient of variation, alias of variation()
        se,          # standard error, alias of sem()
        tab,         # n-way freq table based on FreqTables, NAs are allowed now
        tabprint,    # print a named array output from tab or freqtable
        tabstat,     # compute univariate statistics by subgroups
        substat,     # create a dataarray of univariate statistics of one variable by subgroups
        rowsum,      # compute
        rowstat,     #
        duplicates,  # report, drop, or tag duplicate rows
        destring,    # convert strings to numeric values in a DataArray
        destring!,   # in-place versionof destring
        renvars!,    # change variable names to either lower or upper case
        eform,       # coeftable output to eform
        anova, oneway, # oneway ANOVA
        pwcorr,      # pairwise correlations
        # ranksum,     # ranksum test
        # signrank,    # signed rank test
        ttest,       # t-test
        dir,        # directory listing
        ds,         # filenames according to type, length, or regex
        # getmaxwidth, # maximum length of a string variable
        eltype2,
        nulldeviance,  # nulldeviance for GLM models
        rowpct, colpct, cellpct, chi2, # freqtable functions
        uncategorical!, uncategorize, # functions to reverse CategoricalArray(s) to their original values
        convert_feather, # coonvert Arrow data types in feather file to Julia data types
        identify_condition, identify_condition2 # used to identify conditions in claims files

abstract type Link end
abstract type Formula end
abstract type DictEncoding end


##############################################################################
##
## Load files
##
##############################################################################
include("stella_tools.jl")
include("DataFrame_tools.jl")
include("tab.jl")
include("t-test.jl")
include("anova.jl")
include("show.jl")
include("other.jl")

end
