using DataFrames, Stella, PyCall, JLD

df = load("c:\\Users\\Min\\Dropbox\\Julialang\\mi.jld","data")

desc(df)

countmap(df[:smk])
pool!(df,[:smk])

using DataStructures, GLM

models = OrderedDict()

fm = Vector{Formula}(3)
fm[1] = @formula(mi ~ smk)
fm[2] = @formula(mi ~ smk + sbp)
fm[3] = @formula(mi ~ smk + sbp + ecg)

for i = 1:3
    models[i] = glm(fm[i],df,Bernoulli(),LogitLink())
end

print(coeftable(models[3]))
fieldnames(models[3])
models[3].mf.terms.eterms

depvars = Vector{Symbol}(3)

# dependent variables
for i=1:3
    depvars[i] = models[i].mf.terms.eterms[1]
end

mtitle = nothing
if mtitle == nothing
    mtitle = string.(depvars)
end

coeft = coeftable(models[1])

# collate variables
covariates = Vector{String}()
for i=1:3
    for nm in coeftable(models[i]).rownms
        # println(nm)
        if in(nm, covariates) == false
            println(nm)
            push!(covariates,nm)
        end
    end
end

using PyCall
@pyimport xlsxwriter

wdir = "c:\\data\\Julia"
cd(wdir)

# create a workbook
wb = xlsxwriter.Workbook("summary_logit2.xlsx")

formats = Stella.attach_formats(wb)

t = wb[:add_worksheet]("Summary")

mtitle = ["Model 1: Unadjusted","Model 2: Adjusted for SBP","Model 3: Adjusted for SBP and ECG"]

# determine type of variable
nlev = length(levels(df[:smk]))
otype = "OR"
ci = true
eform = true

num_models = length(models)

r = 0
c = 0

# set column widths
t[:set_column](c,c,40)
t[:set_column](c+1,c+4*(nlev-1),7)

# column title
ctitle = Vector{String}(nlev-1)
ctitle[1] = "smk == 1 (vs smk == 0)"

t[:merge_range](r,c,r+1,c,"Models",formats[:heading])
for i=1:nlev-1
    if ci == true
        t[:merge_range](r,c+1,r,c+4,ctitle[i],formats[:heading])
        t[:merge_range](r+1,c+1,r+1,c+3,string(otype," (95% CI)"),formats[:heading])
        t[:write_string](r+1,c+4,"P-Value",formats[:heading])
    else
        t[:write_string](r,c+1,otype,formats[:heading])
        t[:write_string](r,c+2,"SE",formats[:heading])
        t[:write_string](r,c+3,"Z Value",formats[:heading])
        t[:write_string](r,c+4,"P-Value",formats[:heading])
    end
    c += 4
end

r += 2
c = 0

#println(coeftable(models[3]).rownms)

for i=1:num_models

    # model name
    t[:write_string](r,c,mtitle[i],formats[:model_name])

    # coeftable
    tdata = coeftable(models[i])

    # confint
    tconfint = confint(models[i])

    # estimate index
    k = findfirst(x->x == "smk: 1.0",tdata.rownms)

    # estimates
    for j=1:nlev-1
        # estimates
        if eform == true
            t[:write](r,c+1,exp(tdata.cols[1][k]),formats[:or_fmt])
        else
            t[:write](r,c+1,tdata.cols[1][k],formats[:or_fmt])
        end

        if ci == true

            if eform == true
                # 95% CI Lower
                t[:write](r,c+2,exp(tconfint[k,1]),formats[:cilb_fmt])

                # 95% CI Upper
                t[:write](r,c+3,exp(tconfint[k,2]),formats[:ciub_fmt])
            else
                # 95% CI Lower
                t[:write](r,c+2,tconfint[k,1],formats[:cilb_fmt])

                # 95% CI Upper
                t[:write](r,c+3,tconfint[k,2],formats[:ciub_fmt])
            end
        else
            # SE
            if eform == true
                t[:write](r,c+2,exp(tdata.cols[1][k])*tdata.cols[2][k],formats[:or_fmt])
            else
                t[:write](r,c+2,tdata.cols[1][k],formats[:or_fmt])
            end

            # Z value
            t[:write](r,c+3,tdata.cols[3][k],formats[:or_fmt])

        end

        # P-Value
        t[:write](r,c+4,tdata.cols[4][k].v < 0.001 ? "< 0.001" : tdata.cols[4][k].v ,formats[:p_fmt])

        c += 4

    end

    r += 1
    c = 0
end

wb[:close]()















mglmxls(models,wb,"Models",mtitle=title,eform=true)


# create a worksheet
t = wb[:add_worksheet]("Logit")

# attach formats to the workbook
formats = Stella.attach_formats(wb)

# starting location in the worksheet
r = 0
c = 0

# set column widths
t[:set_column](c,c,40)
t[:set_column](c+1,c+4*3,7)

# headings ---------------------------------------------------------
# if ci == true, Estimate (95% CI) P-Value
# if eform == true, Estimate is OR for logit, IRR for poisson
otype = "OR"
# if eform == true
#     if in(distrib, ["Binomial","Bernoulli"]) && linkfun == "LogitLink"
#         otype = "OR"
#     elseif distrib == "Binomial" && linkfun == "LogLink"
#         otype = "RR"
#     elseif distrib == "Poisson" && linkfun == "LogLink"
#         otype = "IRR"
#     else
#         otype = "exp(Est)"
#     end
# end
t[:merge_range](r,c,r+1,c,"Variable",formats[:heading])
for i=1:3
    t[:merge_range](r,c+1,r,c+4,mtitle[i],formats[:heading])
    t[:merge_range](r+1,c+1,r+1,c+3,string(otype," (95% CI)"),formats[:heading])
    t[:write_string](r+1,c+4,"P-Value",formats[:heading])
    c += 4
end

#------------------------------------------------------------------
r += 2
c = 0
#
# if label_dict != nothing
#     # variable labels
#     varlab = label_dict["variable"]
#
#     # value labels
#     vallab = label_dict["value"]
# end

# go through each variable and construct variable name and value label arrays
tdata = Vector(3)
for i=1:3
    tdata[i] = coeftable(models[i])
end
nrows = length(covariates)
varname = Vector{String}(nrows)
vals = Vector{String}(nrows)
nlev = zeros(Int,nrows)
#nord = zeros(Int,nrow)
for i = 1:nrows
    # variable name
    # parse varname to separate variable name from value
    if contains(covariates[i],":")
        (varname[i],vals[i]) = split(covariates[i],": ")
    else
        varname[i] = covariates[i]
        vals[i] = ""
    end

    # use labels if exist
    # if label_dict != nothing
    #     if vals[i] != ""
    #         valn = vals[i] == "true" ? 1 : parse(Int,vals[i])
    #         if haskey(vallab,varname[i]) && haskey(vallab[varname[i]],valn)
    #             vals[i] = vallab[varname[i]][valn]
    #         end
    #     end
    #     if haskey(varlab,varname[i])
    #         varname[i] = varlab[varname[i]] == "" ? varname[i] : varlab[varname[i]]
    #     end
    # end
end
for i = 1:nrows
    nlev[i] = Stella.countlev(varname[i],varname)
end

eform = true
ci = true

# write table
tconfint = Vector(3)
for i=1:3
    tconfint[i] = confint(models[i])
end
lastvarname = ""

for i = 1:nrows
    # if varname[i] != lastvarname
    #     if nlev[i] > 1
    #         t[:write_string](r,c,varname[i],formats[:heading_left])
    #
    #         if ci == true
    #             t[:write](r,c+1,"",formats[:empty_right])
    #             t[:write](r,c+2,"",formats[:empty_both])
    #             t[:write](r,c+3,"",formats[:empty_left])
    #             t[:write](r,c+4,"",formats[:p_fmt])
    #         else
    #             t[:write](r,c+1,"",formats[:empty_border])
    #             t[:write](r,c+2,"",formats[:empty_border])
    #             t[:write](r,c+3,"",formats[:empty_border])
    #             t[:write](r,c+4,"",formats[:p_fmt])
    #         end
    #         r += 1
    #         t[:write_string](r,c,vals[i],formats[:varname_1indent])
    #
    #     else
    #         if vals[i] != ""
    #             t[:write_string](r,c,string(varname[i]," - ",vals[i]),formats[:heading_left])
    #         else
    #             t[:write_string](r,c,varname[i],formats[:heading_left])
    #         end
    #     end
    # else
    #     t[:write_string](r,c,vals[i],formats[:varname_1indent])
    # end

    # write variable name
    t[:write_string](r,c,varname[i],formats[:model_name])

    for j=1:3

        # find the index number for each coeftable row
        ri = findfirst(x->x == covariates[i],tdata[j].rownms)

        if ri == 0

            # this variable is not in the model
            # print empty cells and then move onto the next model

            t[:write_string](r,c+1,"",formats[:or_fmt])
            t[:write_string](r,c+2,"",formats[:cilb_fmt])
            t[:write_string](r,c+3,"",formats[:ciub_fmt])
            t[:write_string](r,c+4,"",formats[:p_fmt])

            c += 4
            continue

        end

        # estimates
        if eform == true
            t[:write](r,c+1,exp(tdata[j].cols[1][ri]),formats[:or_fmt])
        else
            t[:write](r,c+1,tdata[j].cols[1][ri],formats[:or_fmt])
        end

        if ci == true

            if eform == true
                # 95% CI Lower
                t[:write](r,c+2,exp(tconfint[j][ri,1]),formats[:cilb_fmt])

                # 95% CI Upper
                t[:write](r,c+3,exp(tconfint[j][ri,2]),formats[:ciub_fmt])
            else
                # 95% CI Lower
                t[:write](r,c+2,tconfint[j][ri,1],formats[:cilb_fmt])

                # 95% CI Upper
                t[:write](r,c+3,tconfint[j][ri,2],formats[:ciub_fmt])
            end
        else
            # SE
            if eform == true
                t[:write](r,c+2,exp(tdata[j].cols[1][ri])*tdata[j].cols[2][ri],formats[:or_fmt])
            else
                t[:write](r,c+2,tdata[j].cols[1][ri],formats[:or_fmt])
            end

            # Z value
            t[:write](r,c+3,tdata[j].cols[3][ri],formats[:or_fmt])

        end

        # P-Value
        if varname[i] != lastvarname
            if nlev[i] > 1
                t[:merge_range](r,c+4,r+nlev[i]-1,c+4,tdata[j].cols[4][ri].v < 0.001 ? "< 0.001" : tdata[j].cols[4][ri].v ,formats[:p_fmt])
            else
                t[:write](r,c+4,tdata[j].cols[4][ri].v < 0.001 ? "< 0.001" : tdata[j].cols[4][ri].v ,formats[:p_fmt])
            end
        end

        # update column
        c += 4
    end
    lastvarname = varname[i]
    # update row
    r += 1
    c = 0
end

wb[:close]()
