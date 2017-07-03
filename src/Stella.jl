VERSION >= v"0.5.0" && __precompile__()

module Stella

################################################################################
##
## Dependencies
##
################################################################################

using DataFrames, Distributions, GLM, StatsBase, Glob, Humanize,
        DataStructures, HypothesisTests, NamedArrays, FreqTables

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
        dflist,      # list
        dfmerge,     # merge two dataframes
        dfsample,    # select a sample from a df
        dfappend,    # stack two df's
        univariate,  # univariate statistics in DF
        tabstat,     # unviariate statistics by subgroups
        substat,     # attach univariate stat to original DF by subgroups
        xtile,       # create variable that classify a column into percentiles
        rowsum,      # create variable that sums rows
        strval,      # convert floats or ints to strings
        print,       # print NamedArray
        chisq_2way,  # compute chisquare statistics from na.array
        prepend_spaces, append_spaces, # create fixed length strings
        smallest, largest, # list smallest and largest values in a DA
        toNA,        # set certain values to NA
        pickone,     # create a binary variable that identifies one records in a subgroup
        p5, p10, p25, p50, p75, p90, p95, # percentiles
        cv,          # coefficient of variation, alias of variation()
        se,          # standard error, alias of sem()
        tab,         # n-way freq table based on FreqTables, NAs are allowed now
        tabstat,     # compute univariate statistics by subgroups
        substat,     # create a dataarray of univariate statistics of one variable by subgroups
        rowsum,      # compute
        rowstat,     #
        duplicates,  # report, drop, or tag duplicate rows
        destring,    # convert strings to numeric values in a DataArray
        destring!,   # in-place versionof destring
        renvars!,    # change variable names to either lower or upper case
        recode,      # recode values in a variable
        eform,       # coeftable output to eform
        anova, anovap, # ANOVA table and p-value from ANOVA
        pwcorr,      # pairwise correlations
        ranksum,     # ranksum test
        signrank,    # signed rank test
        ttest,       # t-test
        univariatexls, # output univariate statistics in an excel worksheet
        bivariatexls, # output bivariate statistics in an excel worksheet
        glmxls,     # output GLM output in an excel worksheet
        dfxls,      # output dataframe in an excel file
        dir,        # directory listing
        ds,         # filenames according to type, length, or regex
        getmaxlength # maximum length of a string variable




##############################################################################
##
## Load files
##
##############################################################################

include("stella_tools.jl")
include("DataFrame_tools.jl")
include("Stata_Reader.jl")
include("xlsout.jl")

end
