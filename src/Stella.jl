module Stella

################################################################################
##
## Dependencies
##
################################################################################

using DataFrames, DataAPI, Distributions, StatsBase, StatsAPI, GLM, Survival, 
        NamedArrays, HypothesisTests, DataStructures, FreqTables, 
        LinearAlgebra, Printf, Glob, Dates, Arrow, CategoricalArrays, 
        PooledArrays, PrettyTables, JLD2, Reexport, Plots,
        LogisticROC, TableMetadataTools, SweepOperator, StataFiles

@reexport import StataFiles
@reexport import LogisticROC
@reexport import TableMetadataTools


##############################################################################
##
## Exported methods and types (in addition to everything reexported above)
##
##############################################################################

export  dfcompress, # compress DF
        acompress,  # compress a Vector
        descr,        # list variables with Stata labels and value labels
        ci,          # produce column index from column names
        nmissing,    # returns the number of missing values in an AbstractArray
        lift,        # converts missing values to false in a Boolean array
        filter2,     # returns a boolean vector for use in selecting rows from a DataFrame
        dfsample,    # select a sample from a df
        univariate, univ, # univariate statistics
        tab,         # n-way freq table based on FreqTables, NAs are allowed now
        tabi,        # immediate tab for 2x2 array
        tabstat,     # unviariate statistics by subgroups
        xtile,       # create variable that classify a column into percentiles
        strval,      # convert floats or ints to strings
        chi2test, chi2,    # compute chisquare statistics from na.array
        # prepend_spaces, append_spaces, # create fixed length strings
        smallest, largest, # list smallest and largest values in a DA
        pickone!,     # create a binary variable that identifies one records in a subgroup
        p5, p10, p25, p50, p75, p90, p95, # percentiles
        cv,          # coefficient of variation, alias of variation()
        se,          # standard error, alias of sem()
        tabstat,     # compute univariate statistics by subgroups
        duplicates,  # report, drop, or tag duplicate rows
        destring,    # convert strings to numeric values in a DataArray
        destring!,   # in-place versionof destring
        renvars!,    # change variable names to either lower or upper case
        eform,       # coeftable output to eform
        pwcorr,      # pairwise correlations
        # ranksum,     # ranksum test
        # signrank,    # signed rank test
        ttest,       # t-test
        dir,        # directory listing
        ds,         # filenames according to type, length, or regex
        # getmaxwidth, # maximum length of a string variable
        eltype2,
        rowpct, colpct, cellpct, chi2, # freqtable functions
        categorical!, uncategorical!, uncategorical, # functions to create CategoricalArrays or reverse them to their original values
        # identify_condition, identify_condition2, # used to identify conditions in claims files
        # stats for variables on the same row
        rowtotal, rowfirst, rowlast, rowmean, 
        rowsd, rowmin, rowmax, rowpctile, 
        rowmedian, rowmiss, rownonmiss,
        # subgroup stats
        subtotal!, submean!, subsd!, submin!, submax!, subpctile!,
        submedian!, submiss!, subnonmiss!,
        values!,
        labels!,
        labels2

abstract type Link end
abstract type Formula end

##############################################################################
##
## Load files
##
##############################################################################
include("stella_tools.jl")
include("df_tools.jl")
include("labels.jl")
include("tab.jl")
include("t-test.jl")
include("other.jl")

end
