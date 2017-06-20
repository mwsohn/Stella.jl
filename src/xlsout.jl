using DataFrames, GLM, PyCall, FreqTables

"""
    glmxls(glmout::DataFrames.DataFrameRegressionModel, workbook::PyObject, worksheet::AbstractString; label_dict::Union{Void,Dict}=nothing,eform=false,ci=true, row = 0, col =0)

Outputs a GLM regression table to an excel spreadsheet.
To use this function, `PyCall` is required with a working version python and
a python package called `xlsxwriter` installed. If a label is found for a variable
or a value of a variable in a `label_dict`, the label will be output. Options are:

- `glmout`: returned value from a GLM regression model
- `workbook`: a returned value from xlsxwriter.Workbook() function (see an example below)
- `worksheet`: a string for the worksheet name
- `label_dict`: an option to specify a `label` dictionary (see an example below)
- `eform`: use `eform = true` to get exponentiated estimates, standard errors, or 95% confidence intervals
- `ci`: use `ci = true` (default) to get 95% confidence intervals. `ci = false` will produce standard errors and Z values instead.
- `row`: specify the row of the workbook to start the output table (default = 0 (for row 1))
- `col`: specify the column of the workbook to start the output table (default = 0 (for column A))

# Example 1
This example is useful when one wants to append a worksheet to an existing workbook.
It is responsibility of the user to open a workbook before the function call and close it
to actually create the physical file by close the workbook.

```jldoctest
julia> using PyCall

julia> @pyimport xlsxwriter

julia> wb = xlsxwriter.Workbook("test_workbook.xlsx")
PyObject <xlsxwriter.workbook.Workbook object at 0x000000002A628E80>

julia> glmxls(ols1,wb,"OLS1",label_dict = label)

julia> bivairatexls(df,:incomecat,[:age,:race,:male,:bmicat],wb,"Bivariate",label_dict = label)

Julia> wb[:close]()
```

# Example 2
Alternatively, one can create a spreadsheet file directly. `PyCall` or `@pyimport`
does not need to be called before the function.

```jldoctest
julia> glmxls(ols1,"test_workbook.xlsx","OLS1",label_dict = label)
```

# Example 3
A `label` dictionary is a collection of dictionaries, two of which are `variable` and `value`
dictionaries. The label dictionary can be created as follows:
```jldoctest
julia> label = Dict()
Dict{Any,Any} with 0 entries

julia> label["variable"] = Dict(
    "age" => "Age at baseline",
    "race" => "Race/ethnicity",
    "male" => "Male sex",
    "bmicat" => "Body Mass Index Category")
Dict{String,String} with 4 entries:
  "male"   => "Male sex"
  "race"   => "Race/ethnicity"
  "bmicat" => "Body Mass Index Category"
  "age"    => "Age at baseline"

julia> label["value"] = Dict()
Dict{Any,Any} with 0 entries

julia> label["value"]["race"] = Dict(
    1 => "White",
    2 => "Black",
    3 => "Hispanic",
    4 => "Other")
Dict{Int64,String} with 4 entries:
  4 => "Other"
  2 => "Black"
  3 => "Hispanic"
  1 => "White"

julia> label["value"]["bmicat"] = Dict(
    1 => "< 25 kg/m²",
    2 => "25 - 29.9",
    3 => "20 or over")
Dict{Int64,String} with 3 entries:
  2 => "25 - 29.9"
  3 => "20 or over"
  1 => "< 25 kg/m²"

 julia> label["variable"]["male"]
 "Male sex"

 julia> label["value"]["bmicat"][1]
 "< 25 kg/m²"
 ```
"""
function glmxls(glmout,wbook::PyObject,wsheet::AbstractString;
    label_dict::Union{Void,Dict} = nothing,
    eform::Bool = false, ci = true, row = 0, col = 0)

    modelstr = string(typeof(glmout))
    if !ismatch(r"GLM\.GeneralizedLinearModel",modelstr)
        error("This is not a GLM model output.")
    end

    distrib = replace(modelstr,r".*Distributions\.(Normal|Bernoulli|Binomial|Gamma|Normal|Poisson)\{.*",s"\1")
    linkfun = replace(modelstr,r".*,GLM\.(CauchitLink|CloglogLink|IdentityLink|InverseLink|LogitLink|LogLink|ProbitLink|SqrtLink)\}.*",s"\1")

    # create a worksheet
    t = wb[:add_worksheet](wsheet)

    # attach formats to the workbook
    formats = attach_formats(wbook)

    # starting location in the worksheet
    r = row
    c = col

    # set column widths
    t[:set_column](c,c,40)
    t[:set_column](c+1,c+4,9)

    # headings ---------------------------------------------------------
    # if ci == true, Estimate (95% CI) P-Value
    # if eform == true, Estimate is OR for logit, IRR for poisson
    otype = "Estimate"
    if eform == true
        if distrib == "Binomial" && linkfun == "LogitLink"
            otype = "OR"
        elseif distrib == "Binomial" && linkfun == "LogLink"
            otype = "RR"
        elseif distrib == "Poisson" && linkfun == "LogLink"
            otype = "IRR"
        else
            otype = "exp(Est)"
        end
    end
    t[:write_string](r,c,"Variable",formats[:heading])
    if ci == true
        t[:merge_range](r,c+1,r,c+3,string(otype," (95% CI)"),formats[:heading])
        t[:write_string](r,c+4,"P-Value",formats[:heading])
    else
        t[:write_string](r,c+1,otype,formats[:heading])
        t[:write_string](r,c+2,"SE",formats[:heading])
        t[:write_string](r,c+3,"Z Value",formats[:heading])
        t[:write_string](r,c+4,"P-Value",formats[:heading])
    end

    #------------------------------------------------------------------
    r += 1
    c = col

    if label_dict != nothing
        # variable labels
        varlab = label_dict["variable"]

        # value labels
        vallab = label_dict["value"]
    end

    # go through each variable and construct variable name and value label arrays
    tdata = coeftable(glmout)
    nrows = length(tdata.rownms)
    varname = Vector{String}(nrows)
    vals = Vector{String}(nrows)
    nlev = zeros(Int,nrows)
    #nord = zeros(Int,nrow)
    for i = 1:nrows
    	# variable name
        # parse varname to separate variable name from value
        if contains(tdata.rownms[i],":")
            (varname[i],vals[i]) = split(tdata.rownms[i],": ")
        else
            varname[i] = tdata.rownms[i]
            vals[i] = ""
        end

        # use labels if exist
        if label_dict != nothing
            if vals[i] != ""
                valn = vals[i] == "true" ? 1 : parse(Int,vals[i])
                if haskey(vallab,varname[i]) && haskey(vallab[varname[i]],valn)
                    vals[i] = vallab[varname[i]][valn]
                end
            end
            if haskey(varlab,varname[i])
                varname[i] = varlab[varname[i]] == "" ? varname[i] : varlab[varname[i]]
            end
        end
    end
    for i = 1:nrows
        nlev[i] = countlev(varname[i],varname)
    end

    # write table
    tconfint = confint(glmout)
    lastvarname = ""

    for i = 1:nrows
        if varname[i] != lastvarname
            # output cell boundaries only and go to the next line
            if nlev[i] > 1
                t[:write_string](r,c,varname[i],formats[:heading_left])

                if ci == true
                    t[:write](r,c+1,"",formats[:empty_right])
                    t[:write](r,c+2,"",formats[:empty_both])
                    t[:write](r,c+3,"",formats[:empty_left])
                    t[:write](r,c+4,"",formats[:p_fmt])
                else
                    t[:write](r,c+1,"",formats[:empty_border])
                    t[:write](r,c+2,"",formats[:empty_border])
                    t[:write](r,c+3,"",formats[:empty_border])
                    t[:write](r,c+4,"",formats[:p_fmt])
                end
                r += 1
                t[:write_string](r,c,vals[i],formats[:varname_1indent])

            else
                if vals[i] != ""
                    t[:write_string](r,c,string(varname[i]," - ",vals[i]),formats[:heading_left])
                else
                    t[:write_string](r,c,varname[i],formats[:heading_left])
                end
            end
        else
            t[:write_string](r,c,vals[i],formats[:varname_1indent])
        end

    	# estimates
        if eform == true
    	    t[:write](r,c+1,exp(tdata.cols[1][i]),formats[:or_fmt])
        else
            t[:write](r,c+1,tdata.cols[1][i],formats[:or_fmt])
        end

        if ci == true

            if eform == true
            	# 95% CI Lower
            	t[:write](r,c+2,exp(tconfint[i,1]),formats[:cilb_fmt])

            	# 95% CI Upper
            	t[:write](r,c+3,exp(tconfint[i,2]),formats[:ciub_fmt])
            else
                # 95% CI Lower
            	t[:write](r,c+2,tconfint[i,1],formats[:cilb_fmt])

            	# 95% CI Upper
            	t[:write](r,c+3,tconfint[i,2],formats[:ciub_fmt])
            end
        else
            # SE
            if eform == true
        	    t[:write](r,c+2,exp(tdata.cols[1][i])*tdata.cols[2][i],formats[:or_fmt])
            else
                t[:write](r,c+2,tdata.cols[1][i],formats[:or_fmt])
            end

            # Z value
            t[:write](r,c+3,tdata.cols[3][i],formats[:or_fmt])

        end

        # P-Value
        if varname[i] != lastvarname
            if nlev[i] > 1
                t[:merge_range](r,c+4,r+nlev[i]-1,c+4,tdata.cols[4][i].v < 0.001 ? "< 0.001" : tdata.cols[4][i].v ,formats[:p_fmt])
            else
    	        t[:write](r,c+4,tdata.cols[4][i].v < 0.001 ? "< 0.001" : tdata.cols[4][i].v ,formats[:p_fmt])
            end
        end

        lastvarname = varname[i]

        # update row
        r += 1
    end
end
function glmxls(glmout,wbook::AbstractString,wsheet::AbstractString;
    label_dict::Union{Void,Dict} = nothing,
    eform::Bool = false, ci = true, row = 0, col = 0)

    xlsxwriter = pyimport("xlsxwriter")

    wb = xlsxwriter[:Workbook](wbook)

    glmxls(glmout,wb,wsheet,label_dict=label_dict,eform=eform,ci=ci,row=row,col=col)
end

"""
    bivariatexls(df::DataFrame,colvar::Symbol,rowvars::Vector{Symbol},workbook::PyObject,worksheet::AbstractString; label_dict::Union{Void,Dict}=nothing,row=0,col=0)

 Creates bivariate statistics and appends it in a nice tabular format to an existing workbook.
 To use this function, `PyCall` is required with a working version python and
 a python package called `xlsxwriter` installed.  If a label is found for a variable
 or a value of a variable in a `label_dict`, the label will be output. Options are:

- `df`: a DataFrame
- `colvar`: a categorical variable whose values will be displayed on the columns
- `rowvars`: a Vector of Symbols for variables to be displayed on the rows. Both continuous and categorical variables are allowed.
    For continuous variables, mean and standard deviations will be output and a p-value will be based on an ANOVA test. For categorical variables,
    a r x c table with cell counts and row percentages will be output with a p-value based on a chi-square test.
- `workbook`: a returned value from xlsxwriter.Workbook() function (see an example below)
- `worksheet`: a string for the worksheet name
- `label_dict`: an option to specify a `label` dictionary (see an example below)
- `row`: specify the row of the workbook to start the output table (default = 0 (for row 1))
- `col`: specify the column of the workbook to start the output table (default = 0 (for column A))

# Example 1
This example is useful when one wants to append a worksheet to an existing workbook.
It is responsibility of the user to open a workbook before the function call and close it
to actually create the physical file by close the workbook.

```jldoctest
julia> using PyCall

julia> @pyimport xlsxwriter

julia> wb = xlsxwriter.Workbook("test_workbook.xlsx")
PyObject <xlsxwriter.workbook.Workbook object at 0x000000002A628E80>

julia> glmxls(ols1,wb,"OLS1",label_dict = label)

julia> bivairatexls(df,:incomecat,[:age,:race,:male,:bmicat],wb,"Bivariate",label_dict = label)

Julia> wb[:close]()
```

# Example 2
Alternatively, one can create a spreadsheet file directly. `PyCall` or `@pyimport`
does not need to be called before the function.

```jldoctest
julia> bivariatexls(df,:incomecat,[:age,:race,:male,:bmicat],"test_workbook.xlsx","Bivariate",label_dict = label)
```

# Example 3
A `label` dictionary is a collection of dictionaries, two of which are `variable` and `value`
dictionaries. The label dictionary can be created as follows:

```jldoctest
julia> label = Dict()
Dict{Any,Any} with 0 entries

julia> label["variable"] = Dict(
    "age" => "Age at baseline",
    "race" => "Race/ethnicity",
    "male" => "Male sex",
    "bmicat" => "Body Mass Index Category")
Dict{String,String} with 4 entries:
  "male"   => "Male sex"
  "race"   => "Race/ethnicity"
  "bmicat" => "Body Mass Index Category"
  "age"    => "Age at baseline"

julia> label["value"] = Dict()
Dict{Any,Any} with 0 entries

julia> label["value"]["race"] = Dict(
    1 => "White",
    2 => "Black",
    3 => "Hispanic",
    4 => "Other")
Dict{Int64,String} with 4 entries:
  4 => "Other"
  2 => "Black"
  3 => "Hispanic"
  1 => "White"

julia> label["value"]["bmicat"] = Dict(
    1 => "< 25 kg/m²",
    2 => "25 - 29.9",
    3 => "20 or over")
Dict{Int64,String} with 3 entries:
  2 => "25 - 29.9"
  3 => "20 or over"
  1 => "< 25 kg/m²"

 julia> label["variable"]["male"]
 "Male sex"

 julia> label["value"]["bmicat"][1]
 "< 25 kg/m²"
 ```

"""
function bivariatexls(df::DataFrame,
    colvar::Symbol,
    rowvars::Vector{Symbol},
    wbook::PyObject,
    wsheet::AbstractString;
    label_dict::Union{Void,Dict} = nothing,
    row::Int = 0,
    col::Int = 0)

    if label_dict != nothing
        # variable labels
        varlab = label_dict["variable"]

        # value labels
        vallab = label_dict["value"]
    end

    # create a worksheet
    t = wbook[:add_worksheet](wsheet)

    # attach formats to the workbook
    formats = attach_formats(wbook)

    # starting row and column
    r = row
    c = col

    # drop NAs in colvar
    df2 = df[completecases(df[[colvar]]),:]

    # number of columns
    # column values
    collev = freqtable(df2,colvar)
    nlev = length(collev.array)
    colnms = names(collev,1)
    coltot = sum(collev.array,1)

    # set column widths
    t[:set_column](c,c,40)
    t[:set_column](c+1,c+(nlev+1)*2+1,9)

    # create heading
    # column variable name
    # It uses three rows
    t[:merge_range](r,c,r+2,c,"Variable",formats[:heading])

    # 1st row = variable name
    colvname = string(colvar)
    if label_dict != nothing
        if haskey(varlab,colvname)
            colvname = varlab[colvname]
        end
    end
    t[:merge_range](r,c+1,r,c+(nlev+1)*2+1,colvname,formats[:heading])

    # 2nd and 3rd rows
    r += 1

    t[:merge_range](r,1,r,2,"All",formats[:heading])
    t[:write_string](r+1,1,"N",formats[:n_fmt_right])
    t[:write_string](r+1,2,"(%)",formats[:pct_fmt_parens])

    c = 3
    for i = 1:nlev

        # value label
        vals = string(colnms[i])
        if label_dict != nothing
            if haskey(vallab,colnms[i])
                vals = vallab[colnms[i]]
            end
        end

        t[:merge_range](r,c+(i-1)*2,r,c+(i-1)*2+1,vals,formats[:heading])
        t[:write_string](r+1,c+(i-1)*2,"N",formats[:n_fmt_right])
        t[:write_string](r+1,c+(i-1)*2+1,"(%)",formats[:pct_fmt_parens])
    end

    # P-value
    t[:merge_range](r,c+nlev*2,r+1,c+nlev*2,"P-Value",formats[:heading])

    # total
    c = 0
    r += 2
    t[:write_string](r,c,"All",formats[:heading_left])
    x = freqtable(df2,colvar)
    tot = sum(x.array)
    t[:write](r,c+1,tot,formats[:n_fmt_right])
    t[:write](r,c+2,1.0,formats[:pct_fmt_parens])
    for i = 1:nlev
        t[:write](r,c+i*2+1,x.array[i],formats[:n_fmt_right])
        t[:write](r,c+i*2+2,x.array[i]/tot,formats[:pct_fmt_parens])
    end
    t[:write](r,c+(nlev+1)*2+1,"",formats[:empty_border])

    # covariates
    c = 0
    r += 1
    for varname in rowvars

        # print the variable name
        vars = string(varname)
        if label_dict != nothing
            if haskey(varlab,vars)
                vars = varlab[vars]
            end
        end

        # determine if varname is categorical or continuous
        if typeof(df2[varname]) <: PooledDataArray || eltype(df2[varname]) == String
            # categorial
            df3=df2[completecases(df2[[varname]]),[varname,colvar]]
            x = freqtable(df3,varname,colvar)
            rowval = names(x,1)
            rowtot = sum(x.array,2)

            # variable name
            # if there only two levels and one of the values is 1 or true
            # just output the frequency and percentage of the 1/true row
            t[:write_string](r,c,vars,formats[:heading_left])

            for i = 1:nlev+1
                t[:write_string](r,c+(i-1)*2+1,"",formats[:empty_right])
                t[:write_string](r,c+(i-1)*2+2,"",formats[:empty_left])
            end
            t[:write_string](r,c+(nlev+1)*2+1,"",formats[:empty_border])

            r += 1
            for i = 1:length(rowval)
                # row value
                vals = string(rowval[i])
                if label_dict != nothing
                    vars = string(varname)
                    if haskey(vallab, vars) && haskey(vallab[vars],rowval[i])
                        vals = vallab[vars][rowval[i]]
                    end
                end
                t[:write_string](r,c,vals,formats[:varname_1indent])

                # row total
                t[:write](r,c+1,rowtot[i],formats[:n_fmt_right])
                t[:write](r,c+2,rowtot[i]/tot,formats[:pct_fmt_parens])

                for j = 1:nlev
                    t[:write](r,c+j*2+1,x.array[i,j],formats[:n_fmt_right])
                    t[:write](r,c+j*2+2,x.array[i,j]/rowtot[i],formats[:pct_fmt_parens])
                end

                # p-value - output only once
                if length(rowval) == 1
                    t[:write](r,c+(nlev+1)*2,chisq_2way(x)[3],formats[:p_fmt])
                elseif i == 1
                    pval = chisq_2way(x)[3]
                    t[:merge_range](r,c+(nlev+1)*2+1,r+length(rowval)-1,c+(nlev+1)*2+1,pval < 0.001 ? "< 0.001" : pval,formats[:p_fmt])
                end
                r += 1
            end
        else
            # continuous variable
            df3=df2[completecases(df2[[varname]]),[varname,colvar]]
            y = tabstat(df3,varname,colvar)

            # variable name
            t[:write_string](r,c,string(vars,", mean (SD)"),formats[:heading_left])

            # All
            t[:write](r,c+1,mean(df3[varname]),formats[:f_fmt_right])
            t[:write](r,c+2,std(df3[varname]),formats[:f_fmt_left_parens])

            # colvar levels
            for i = 1:nlev
                if i <= size(y,1) && y[i,:N] > 1
                    t[:write](r,c+i*2+1,y[i,:mean],formats[:f_fmt_right])
                    t[:write](r,c+i*2+2,y[i,:sd],formats[:f_fmt_left_parens])
                else
                    t[:write](r,c+i*2+1,"",formats[:f_fmt_right])
                    t[:write](r,c+i*2+2,"",formats[:f_fmt_left_parens])
                end
            end
            if size(y,1) > 1
                pval = anovap(df3,varname,colvar)
                t[:write](r,c+(nlev+1)*2+1,pval < 0.001 ? "< 0.001" : pval,formats[:p_fmt])
            else
                t[:write](r,c+(nlev+1)*2+1,"",formats[:p_fmt])
            end

            r += 1
        end
    end
end
function bivariatexls(df::DataFrame,
    colvar::Symbol,
    rowvars::Vector{Symbol},
    wbook::AbstractString,
    wsheet::AbstractString;
    label_dict::Union{Void,Dict} = nothing,
    row::Int = 0,
    col::Int = 0)

    xlsxwriter = pyimport("xlsxwriter")

    wb = xlsxwriter[:Workbook](wbook)

    bivariatexls(df,colvar,rowvars,wb,wsheet,label_dict=label_dict,row=row,col=col)
end


"""
    univariatexls(df::DataFrame,contvars::Vector{Symbol},workbook::PyObject,worksheet::AbstractString; label_dict::Union{Void,Dict}=nothing,row=0,col=0)

Creates univariate statistics for a vector of continuous variable and
appends it to an existing workbook.
To use this function, `PyCall` is required with a working version python and
a python package called `xlsxwriter` installed.  If a label is found for a variable
in a `label_dict`, the label will be output. Options are:

- `df`: a DataFrame
- `contvars`: a vector of continuous variables
- `workbook`: a returned value from xlsxwriter.Workbook() function (see an example below)
- `worksheet`: a string for the worksheet name
- `label_dict`: an option to specify a `label` dictionary (see an example below)
- `row`: specify the row of the workbook to start the output table (default = 0 (for row 1))
- `col`: specify the column of the workbook to start the output table (default = 0 (for column A))

# Example 1
This example is useful when one wants to append a worksheet to an existing workbook.
It is responsibility of the user to open a workbook before the function call and close it
to actually create the physical file by close the workbook.

```jldoctest
julia> using PyCall

julia> @pyimport xlsxwriter

julia> wb = xlsxwriter.Workbook("test_workbook.xlsx")
PyObject <xlsxwriter.workbook.Workbook object at 0x000000002A628E80>

julia> glmxls(ols1,wb,"OLS1",label_dict = label)

julia> bivairatexls(df,:incomecat,[:age,:race,:male,:bmicat],wb,"Bivariate",label_dict = label)

julia> univariatexls(df,[:age,:income_amt,:bmi],wb,"Univariate",label_dict = label)

Julia> wb[:close]()
```

# Example 2
Alternatively, one can create a spreadsheet file directly. `PyCall` or `@pyimport`
does not need to be called before the function.

```jldoctest
julia> univariatexls(df,[:age,:income_amt,:bmi],"test_workbook.xlsx","Bivariate",label_dict = label)
```

# Example 3
A `label` dictionary is a collection of dictionaries, two of which are `variable` and `value`
dictionaries. The label dictionary can be created as follows:

```jldoctest
julia> label = Dict()
Dict{Any,Any} with 0 entries

julia> label["variable"] = Dict(
    "age" => "Age at baseline",
    "race" => "Race/ethnicity",
    "male" => "Male sex",
    "bmicat" => "Body Mass Index Category")
Dict{String,String} with 4 entries:
  "male"   => "Male sex"
  "race"   => "Race/ethnicity"
  "bmicat" => "Body Mass Index Category"
  "age"    => "Age at baseline"

 julia> label["variable"]["male"]
 "Male sex"
 ```

"""
function univariatexls(df::DataFrame,contvars::Vector{Symbol},wbook::PyObject,wsheet::AbstractString;
    label_dict::Union{Void,Dict}=nothing, row = 0, col = 0)

    if label_dict != nothing
        # variable labels
        varlab = label_dict["variable"]
    end

    # create a worksheet
    t = wbook[:add_worksheet](wsheet)

    # attach formats to the workbook
    formats = attach_formats(wbook)

    # starting row and column
    r = row
    c = col

    # column width
    t[:set_column](0,0,20)
    t[:set_column](1,length(contvars),12)

    # output the row names
    t[:write_string](r,c,"Statistic",formats[:heading])
    t[:write_string](r+1,c,"N Total",formats[:heading_left])
    t[:write_string](r+2,c,"N Miss",formats[:heading_left])
    t[:write_string](r+3,c,"N Used",formats[:heading_left])
    t[:write_string](r+4,c,"Sum",formats[:heading_left])
    t[:write_string](r+5,c,"Mean",formats[:heading_left])
    t[:write_string](r+6,c,"SD",formats[:heading_left])
    t[:write_string](r+7,c,"Variance",formats[:heading_left])
    t[:write_string](r+8,c,"Minimum",formats[:heading_left])
    t[:write_string](r+9,c,"P25",formats[:heading_left])
    t[:write_string](r+10,c,"Median",formats[:heading_left])
    t[:write_string](r+11,c,"P75",formats[:heading_left])
    t[:write_string](r+12,c,"Maximum",formats[:heading_left])
    t[:write_string](r+13,c,"Skewness",formats[:heading_left])
    t[:write_string](r+14,c,"Kurtosis",formats[:heading_left])
    t[:write_string](r+15,c,"Smallest",formats[:heading_left])
    t[:write_string](r+16,c,"",formats[:heading_left])
    t[:write_string](r+17,c,"",formats[:heading_left])
    t[:write_string](r+18,c,"",formats[:heading_left])
    t[:write_string](r+19,c,"",formats[:heading_left])
    t[:write_string](r+20,c,"Largest",formats[:heading_left])
    t[:write_string](r+21,c,"",formats[:heading_left])
    t[:write_string](r+22,c,"",formats[:heading_left])
    t[:write_string](r+23,c,"",formats[:heading_left])
    t[:write_string](r+24,c,"",formats[:heading_left])

    col = 1
    for vsym in contvars

        # if there is a label dictionary, pick up the variable label
        varstr = string(vsym)
        if label_dict != nothing && haskey(varlab,varstr)
            varstr = varlab[varstr] == "" ? varstr : varlab[varstr]
        end

        t[:write_string](0,col,varstr,formats[:heading])
        u = univariate(df,vsym)
        for j = 1:14
            if j<4
                fmttype = :n_fmt
            else
                fmttype = :p_fmt
            end
            t[:write](j,col,u[j,:Value],formats[fmttype])
        end
        smallest=Stella.smallest(df[vsym])
        if eltype(df) <: Integer
            fmttype = :n_fmt
        else
            fmttype = :p_fmt
        end
        for j = 1:5
            t[:write](j+14,col,smallest[j],formats[fmttype])
        end
        largest=Stella.largest(df[vsym])
        for j = 1:5
            t[:write](j+19,col,largest[j],formats[fmttype])
        end
        col += 1
    end
end
function univariatexls(df::DataFrame,contvars::Vector{Symbol},wbook::AbstractString,wsheet::AbstractString;
    label_dict::Union{Void,Dict}=nothing, row = 0, col = 0)

    xlsxwriter=pyimport("xlsxwriter")

    wb = wlsxwriter[:Workbook](wbook)

    univariatexls(df,contvars,wb,wsheet,label_dict=label_dict,row=row,col=col)
end


"""
    dfxls(df::DataFrame,workbook::PyObject; worksheet::AbstractString = "Data1",nrow = 500, start = 1, row=0,col=0)

 To use this function, `PyCall` is required with a working version python and
 a python package called `xlsxwriter` installed. Options are:

- `df`: a DataFrame
- `workbook`: a returned value from xlsxwriter.Workbook() function (see an example below)
- `worksheet`: a string for the worksheet name (default: "Data1")
- `start`: specify the row number from which data will be output (default: 1)
- `nrows`: specify the number of rows to output (default: 500). If nows = 0, the entire dataframe will be output.
- `row`: specify the row of the workbook to start the output table (default = 0 (for row 1))
- `col`: specify the column of the workbook to start the output table (default = 0 (for column A))

# Example 1
This example is useful when one wants to append a worksheet to an existing workbook.
It is responsibility of the user to open a workbook before the function call and close it
to actually create the physical file by close the workbook.

```jldoctest
julia> using PyCall

julia> @pyimport xlsxwriter

julia> wb = xlsxwriter.Workbook("test_workbook.xlsx")
PyObject <xlsxwriter.workbook.Workbook object at 0x000000002A628E80>

julia> glmxls(ols1,wb,"OLS1",label_dict = label)

julia> bivairatexls(df,:incomecat,[:age,:race,:male,:bmicat],wb,"Bivariate",label_dict = label)

julia> univariatexls(df,[:age,:income_amt,:bmi],wb,"Univariate",label_dict = label)

julia> dfxls(df,wb,"dataframe",nrows = 0)

Julia> wb[:close]()
```

# Example 2
Alternatively, one can create a spreadsheet file directly. `PyCall` or `@pyimport`
does not need to be called before the function.

```jldoctest
julia> dfxls(df,"test_workbook.xlsx","df",nrows = 0)
```

"""
function dfxls(df::DataFrame,
    wbook::PyObject;
    worksheet::AbstractString = "Data1",
    rows::Int64 = 500, start::Int64 = 1, col::Int64 = 0, row::Int64 = 0)

    # create a worksheet
    t = wbook[:add_worksheet](worksheet)

    # attach formats to the workbook
    formats = attach_formats(wbook)

    # starting row and column
    c = col

    typ = eltypes(df)
    varnames = names(df)

    # if rows = 0, output the full data set
    if nrows == 0
        start = 1
        nrows = size(df,1)
    else
        nrows = (start+nrows-1) < size(df,1) ? brows : (size(df,1)-start+1)
    end

    for i in 1:size(df,2)

        r = row
        t[:set_column](c,c,12)
        t[:write_string](r,c,string(varnames[i]),formats[:heading])
        r += 1

        for j in start:(start+nrows-1)

            if isna(df[j,i]) && typ[i] <: Number
                t[:write_string](r,c,"NA",formats[:n_fmt])
            elseif isna(df[j,i]) && typ[i] <: AbstractString
                t[:write_string](r,c,"",formats[:text])
            elseif typ[i] <: Integer
                t[:write](r,c,df[j,i],formats[:n_fmt])
            elseif typ[i] <: AbstractFloat
                t[:write](r,c,df[j,i],formats[:f_fmt])
            else
                t[:write](r,c,df[j,i],formats[:text])
            end

            r += 1
        end
        c += 1
    end

    wb[:close]()
end
function dfxls(df::DataFrame,
    wbook::AbstractString;
    worksheet::AbstractString = "Data1",
    nrows::Int64 = 500, start::Int64 = 1, col::Int64 = 0, row::Int64 = 0)

    # import xlsxwriter
    XlsxWriter = pyimport("xlsxwriter")

    # create a workbook
    wb = XlsxWriter[:Workbook](wbook)

    dfxls(df,wb,worksheet = worksheet, nrows = nrows, start = start, col = col, row = row)
end

function newfilename(filen::AbstractString)
    while (isfile(filen))
        # separate the name into three parts, basename (number).extension
        (basename,ext) = splitext(filen)
        m = match(r"(.*)(\(.*\))$",basename)
        if m == nothing
            filen = string(filen," (1)")
        else
            m2 = parse(Int,replace(m[2],r"\((.*)\)",s"\1"))
            filen = string(m[1]," (",m2+1,")")
        end
    end
    return filen
end

function countlev(str::AbstractString,sarray::Vector{String})
    k = 0
    for i = 1:length(sarray)
        if str == sarray[i]
            k += 1
        end
    end
    return k
end


format_defs = Dict()

format_defs[:heading] = Dict(
	"bold" => true,
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "center",
	"border" => true
)

format_defs[:text] = Dict(
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "left",
	"border" => true
)


format_defs[:heading_right] = Dict(
	"bold" => true,
	"font" => "Arial",
	"size" => 10,
	"valign" => "vcenter",
	"align" => "right",
	"left" => true, "bottom" => true, "top" => true
)

format_defs[:heading_left] = Dict(
	"bold" => true,
	"font" => "Arial",
	"size" => 10,
	"valign" => "vcenter",
	"align" => "left",
	"right" => true, "bottom" => true, "top" => true
)

format_defs[:model_name] = Dict(
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "left",
	"border" => true
)

format_defs[:varname_1indent] = Dict(
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "left",
	"border" => true,
	"indent" => 1
)


format_defs[:varname_2indent] = Dict(
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "left",
	"border" => true,
	"indent" => 2
)


format_defs[:n_fmt_right] = Dict(
	"num_format" => "#,##0",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "right",
	"left" => true, "bottom" => true, "top" => true
)


format_defs[:n_fmt_left_parens] = Dict(
	"num_format" => "(#,##0)",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "right",
	"right" => true, "bottom" => true, "top" => true
)

format_defs[:f_fmt_right] = Dict(
	"num_format" => "#,##0.00",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "right",
	"left" => true, "bottom" => true, "top" => true
)

format_defs[:f_fmt] = Dict(
	"num_format" => "#,##0.00",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "right",
	"border" => true
)


format_defs[:n_fmt] = Dict(
	"num_format" => "#,##0",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "right",
	"border" => true
)


format_defs[:f_fmt_left_parens] = Dict(
	"num_format" => "(#,##0.00)",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "left",
	"right" => true, "bottom" => true, "top" => true
)


format_defs[:pct_fmt_parens] = Dict(
	"num_format" => "(0.00%)",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "left",
	"right" => true, "bottom" => true, "top" => true
)


format_defs[:pct_fmt] = Dict(
	"num_format" => "0.00%",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "right",
	"border" => true
)

format_defs[:or_fmt] = Dict(
	"num_format" => "0.000",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "right",
	"left" => true, "bottom" => true, "top" => true
)

format_defs[:cilb_fmt] = Dict(
	"num_format" => "(0.000 -",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "right",
	"bottom" => true, "top" => true
)

format_defs[:ciub_fmt] = Dict(
	"num_format" => "0.000)",
	"font" => "Arial",
	"font_size" => 9,
	"valign" => "vcenter",
	"align" => "left",
	"right" => true, "bottom" => true, "top" => true
)


format_defs[:or_fmt_red] = Dict(
	"num_format" => "0.000",
	"font" => "Arial",
	"font_size" => 9,
	"font_color" => "red",
	"valign" => "vcenter",
	"align" => "right",
	"left" => true, "bottom" => true, "top" => true
)

format_defs[:cilb_fmt_red] = Dict(
	"num_format" => "(0.000 -",
	"font" => "Arial",
	"font_size" => 9,
	"font_color" => "red",
	"valign" => "vcenter",
	"align" => "right",
	"bottom" => true, "top" => true
)

format_defs[:ciub_fmt_red] = Dict(
	"num_format" => "0.000)",
	"font" => "Arial",
	"font_size" => 9,
	"font_color" => "red",
	"valign" => "vcenter",
	"align" => "left",
	"right" => true, "bottom" => true, "top" => true
)


format_defs[:or_fmt_bold] = Dict(
	"num_format" => "0.000",
	"font" => "Arial",
	"font_size" => 9,
	"bold" => true,
	"valign" => "vcenter",
	"align" => "right",
	"left" => true, "bottom" => true, "top" => true
)

format_defs[:cilb_fmt_bold] = Dict(
	"num_format" => "(0.000 -",
	"font" => "Arial",
	"font_size" => 9,
	"bold" => true,
	"valign" => "vcenter",
	"align" => "right",
	"bottom" => true, "top" => true
)

format_defs[:ciub_fmt_bold] = Dict(
	"num_format" => "0.000)",
	"font" => "Arial",
	"font_size" => 9,
	"bold" => true,
	"valign" => "vcenter",
	"align" => "left",
	"right" => true, "bottom" => true, "top" => true
)

format_defs[:p_fmt] = Dict(
	"num_format" => "0.000",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "right",
	"border" => true
)

format_defs[:p_fmt2] = Dict(
	#"num_format" => "0.000",
	"font" => "Arial",
	"size" => 9,
	"valign" => "vcenter",
	"align" => "right",
	"border" => true
)


format_defs[:p_fmt_red] = Dict(
	"num_format" => "0.000",
	"font" => "Arial",
	"font_size" => 9,
	"font_color" => "red",
	"valign" => "vcenter",
	"align" => "right",
	"border" => true
)

format_defs[:p_fmt2_red] = Dict(
	#"num_format" => "0.000",
	"font" => "Arial",
	"font_size" => 9,
	"font_color" => "red",
	"valign" => "vcenter",
	"align" => "right",
	"border" => true
)


format_defs[:p_fmt_bold] = Dict(
	"num_format" => "0.000",
	"font" => "Arial",
	"font_size" => 9,
	"bold" => true,
	"valign" => "vcenter",
	"align" => "right",
	"border" => true
)

format_defs[:p_fmt2_bold] = Dict(
	#"num_format" => "0.000",
	"font" => "Arial",
	"font_size" => 9,
	"bold" => true,
	"valign" => "vcenter",
	"align" => "right",
	"border" => true
)

format_defs[:empty_border] = Dict(
	#"num_format" => "0.000",
	"valign" => "vcenter",
	"border" => true
)

format_defs[:empty_left] = Dict(
	#"num_format" => "0.000",
	"valign" => "vcenter",
	"right" => true, "bottom" => true, "top" => true
)

format_defs[:empty_right] = Dict(
	#"num_format" => "0.000",
	"valign" => "vcenter",
	"left" => true, "bottom" => true, "top" => true
)


format_defs[:empty_both] = Dict(
	#"num_format" => "0.000",
	"valign" => "vcenter",
	"bottom" => true, "top" => true
)


function attach_formats(workbook;formats::Dict = format_defs)
	newfmts = Dict()
	for key in keys(formats)
		newfmts[key] = workbook[:add_format](formats[key])
	end
	return newfmts
end
