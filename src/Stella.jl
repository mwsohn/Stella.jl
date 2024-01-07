module Stella

################################################################################
##
## Dependencies
##
################################################################################

using DataFrames, DataAPI, Distributions, StatsBase, StatsAPI, GLM, NamedArrays, HypothesisTests,
        DataStructures, FreqTables, ReadStat, Labels, LinearAlgebra, Printf, 
        Glob, Dates, Arrow, CategoricalArrays, PooledArrays, PrettyTables, JLD2

##############################################################################
##
## Exported methods and types (in addition to everything reexported above)
##
##############################################################################

export  read_stata,  # read stata 13 and 14 files into DF
        read_stata!, # read stata 13 and 14 files into DF
        dfcompress, # compress DF
        acompress,  # compress a Vector
        desc,        # list variables with Stata labels and value labels
        ci,          # produce column index from column names
        nmissing,    # returns the number of missing values in an AbstractArray
        lift,        # converts missing values to false in a Boolean array
        eq, lt, le, gt, ge, # alternative opreator functions that takes care of missing values
        dflist,      # list
        dfmerge,     # merge two dataframes
        dfsample,    # select a sample from a df
        univariate, univ, # univariate statistics in DF
        tabstat,     # unviariate statistics by subgroups
        substat,     # attach univariate stat to original DF by subgroups
        xtile,       # create variable that classify a column into percentiles
        strval,      # convert floats or ints to strings
        chi2test, chi2,    # compute chisquare statistics from na.array
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
        # nulldeviance,  # nulldeviance for GLM models
        rowpct, colpct, cellpct, chi2, # freqtable functions
        categorical!, uncategorical!, uncategorize, # functions to create CategoricalArrays or reverse them to their original values
        identify_condition, identify_condition2 # used to identify conditions in claims files

abstract type Link end
abstract type Formula end

##############################################################################
##
## Load files
##
##############################################################################
include("alt_operators.jl")
include("stella_tools.jl")
include("DataFrame_tools.jl")
#include("labels.jl")
include("tab.jl")
include("t-test.jl")
include("anova.jl")
#include("show.jl")
include("other.jl")

end
