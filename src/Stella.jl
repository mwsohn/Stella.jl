module Stella

################################################################################
##
## Dependencies
##
################################################################################

using DataFrames, DataAPI, Distributions, StatsBase, StatsAPI, GLM, Survival, 
        NamedArrays, HypothesisTests, DataStructures, FreqTables, 
        LinearAlgebra, Printf, Glob, Dates, Arrow, CategoricalArrays, 
        PooledArrays, PrettyTables, JLD2, Reexport, Plots

@reexport using StataFiles
@reexport using LogisticROC
@reexport using TableMetadataTools
@reexport using AnalysisOfVariance
@reexport using FreqTools


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
        tabstat,     # unviariate statistics by subgroups
        xtile,       # create variable that classify a column into percentiles
        strval,      # convert floats or ints to strings
        chi2test,    # compute chisquare statistics from na.array
        # prepend_spaces, append_spaces, # create fixed length strings
        smallest, largest, # list smallest and largest values in a DA
        pickone,     # create a binary variable that identifies one record in a subgroup
        keepfirst, keeplast,     # creates a DF that contain one record from a subgroup
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
        labels2,
        sasdate, # convert Julia Date to SAS date
        juliadate, # convert SAS date to Julia Date
        st2ncc

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
include("t-test.jl")
include("other.jl")

end
