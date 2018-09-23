# Stella

A package to provide functions that mimic some basic Stata functionality.

## Installation

```
  ] add https://github.com/mwsohn/Stella.jl
```

## List of Functions

Here's the list of functions. Documentation for individual functions can be accessed using the Julia help (e.g., ?read_stata at the Julia REPL prompt).

### DataFrame tools

- read_stata - This is a wrapper for `read_dta` function in ReadStat package that produces a dataframe as output.
- desc - Prints the list of variables, similar to `showcols`. It prints variable and value labels with `labels` option.
- dfcompress! - Compress all numeric variables in a dataframe by converting them to the type that will occupy the least memory.
- acompress - Compress one abstract array
- duplicates - Identifies, deletes, or tag duplicates in a dataframe
- pickone - Identifies only one record in repeated data similar to `egen tag` function in stata
- dfsample - Creates a dataframe with a randomly selected sample from the input dataframe
- dfmerge - Merges two dataframes and creates a `:_merge` variable that indicates data sources, similar to Stata `merge` command
- classify - Creates a classification variable with the specified cutoff values
- xtile - Creates a ntile classification variable, similar to Stata's `xtile` command
- destring - Converts a string vector to a numeric vector
- rowstat - Creates a vector of statistics from multiple rows in a dataframe, analogous to Stata's `egen` row commands
- substat - Creates a vector of statistics for a subgroup of rows, analogous to Stata's `egen, by` commands

### Statitical utilities

- pwcorr - Produces pairwise correlations
- anova - One-way analysis of variance
- univariate - Univariate statistics, output to a DataFrame
- ttest - T-test statistics
- tabstat - Creates statistics by categories using functions provided, similar to Stata's tabstat command
- eform - Produces exponentiated GLM regression output such as coefficients and standard errors
