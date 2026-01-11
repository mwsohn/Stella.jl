# Stella

A package of functions that mimic some basic Stata functionality.

## Installation

```
  ] add https://github.com/mwsohn/Stella.jl
```

## List of Functions

Here's the list of functions. Documentation for individual functions can be accessed using the Julia help (e.g., ?read_stata at the Julia REPL prompt).

### DataFrame tools

- acompress - Compress a vector
- classify - Creates a classification variable with the specified cutoff values
- descr - Prints the list of variables, similar to `showcols` 
- dfcompress - Compress all numeric variables in a dataframe by converting them to the type that will occupy the least memory.
- dfmerge - Merges two dataframes and creates a `:_merge` variable that indicates data sources, similar to Stata `merge` command
- dfsample - Creates a dataframe with a randomly selected sample from the input dataframe
- duplicates - Identifies, deletes, or tag duplicates in a dataframe
- keepfirst - Keep the first `n` rows from the same group in grouped data (default: `n = 1`)
- keeplast - Keep the last `n` rows from the same group in grouped data (default: `n = 1`)
- pickone - Identifies only one record in repeated data similar to `egen tag` function in Stata
- renvars! - Change case for all variables (default: )
- xtile - Creates a ntile classification variable, similar to Stata's `xtile` command

### Row functions that work on a column of values in a subgroup and attaches return values to the original DataFrame

- rowfirst - Returns the first non-missing value
- rowlast - Returns the last non-missing value
- rowmax - Returns the maximum non-missing value
- rowmin - Returns the minimum non-missing value
- rowmean - Returns the mean of non-missing values
- rowmedian - Returns the median of non-missing values
- rowmiss - Returns the number of missing values
- rownonmiss - Returns the number of non-missing values
- rowsd - Returns the standard deviation of non-missing values
- rowpctile - Returns the percentile value of non-missing values (default: `p = 0.5`)
- rowtotal - Returns the sum of non-missing values

### SubDataFrame functions that work on a subgroup vector in a grouped data. A variable is automatically created with "_total", "_mean", "_sd", etc, appended after the variable to generate statistics for.

- subfirst! - Returns the first non-missing value
- sublast! - Returns the last non-missing value
- submax! - Returns the maximum non-missing value
- submin! - Returns the minimum non-missing value
- submean! - Returns the mean of non-missing values
- submedian! - Returns the median of non-missing values
- submiss! - Returns the number of missing values
- subnonmiss! - Returns the number of non-missing values
- subsd! - Returns the standard deviation of non-missing values
- subpctile! - Returns the percentile value of non-missing values (default: `p = 0.5`)
- subtotal! - Returns the sum of non-missing values

### Statitical utilities for interactive use

- pwcorr - Produces pairwise correlations
- univ - Univariate statistics, output to a DataFrame
- ttest - T-test statistics
- tabstat - Creates statistics by categories using functions provided, similar to Stata's tabstat command
- eform - Produces exponentiated GLM regression output such as coefficients and standard errors
