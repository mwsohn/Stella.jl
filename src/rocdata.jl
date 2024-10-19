# adapted from _roc function in ROC.jl to speed up the computation

function rocdata(scores, labels)

    # function that returns number of trues and falses
    sumn(v) = (s = sum(v); return (s, length(v) - s))

    # thresholds are sorted list of unique values in scores
    thresholds = vcat(1.0,sort(unique(scores), rev=true))

    # turn labels into a boolean vector
    labels = labels .== 1.0

    # total number of positives (trues)
    P = sum(labels)

    # number of negatives (falses)
    N = length(labels) - P

    # number of thresholds
    n_thresholds = length(thresholds)

    # allocate memory for ROCData vectors
    TP = Array{Int}(undef, n_thresholds) # true positives
    TN = Array{Int}(undef, n_thresholds) # true negatives
    FP = Array{Int}(undef, n_thresholds) # false positives
    FN = Array{Int}(undef, n_thresholds) # false negatives
    FPR = Array{Float64}(undef, n_thresholds)
    TPR = Array{Float64}(undef, n_thresholds)

    for (i, threshold) in enumerate(thresholds)
        mask = scores .>= threshold
        # predicted_positive = labels[mask]
        # predicted_negative = labels[.!mask]
        # TP[i] = sum(predicted_positive)
        # TN[i] = sum(.!predicted_negative)
        # FP[i] = length(predicted_positive) - TP[i]
        # FN[i] = length(predicted_negative) - TN[i]
        (TP[i], FP[i]) = sumnn(labels[mask])
        (TN[i], FN[i]) = sumnn(labels[.!mask])
        FPR[i] = FP[i] / (FP[i] + TN[i])
        TPR[i] = TP[i] / (TP[i] + FN[i])
    end
    ROCData{eltype(thresholds)}(thresholds, P, N, TP, TN, FP, FN, FPR, TPR)
end

auc(scores, labels) = ROC.AUC(rocdata(scores,labels))
function auc(rr::ROCData)
    println("Number of observations = ", rr.P + rr.N)
    println("Area under ROC curve   = ", AUC(rr))
end

import Plots.plot
function plot(rocdata::ROCData)
    plot(rocdata.FPR, rocdata.TPR, linetype=:steppre, legend=false, xlabel="1 - specificity", ylabel="sensitivity")
    plot!(collect(0:1:1), collect(0:1:1), legend=false, lc = :black)
end


# function logrank(df, event, by)

#     # number of groups
#     # ba = completecases(df[!,[event,by]])
#     # df2 = df[ba,:]
#     groups = sort(unique(df[!,by]))
#     n_groups = length(groups)

#     # perform Kaplan-Meier analysis
#     km = Vector{KaplanMeier{Float64, Int64}}(undef,n_groups)
#     times = zeros(Int64,n_groups)

#     for (i,v) in enumerate(groups)
#         km[i] = fit(KaplanMeier, df[ df[!,by] .== v, event])
#         times[i] = km[i].events.time[end]
#     end

#     # lengths
#     ntimes = maximum(times)

#     # events
#     events = zeros(Int64, n_groups, ntimes)
#     for i in 1:n_groups
#         events[i, km[i].events.time] .= km[i].events.nevents
#     end

#     # N at risk
#     atrisk = zeros(Int64, n_groups, ntimes)
#     for i in 1:n_groups
#         atrisk[i, 1] = km[i].events.natrisk[1]
        
#         for j in 2:ntimes
#             atrisk[i, j] = atrisk[i, j-1] - events[i, j-1]
#         end
#     end

#     # Observed events
#     O = vec(sum(events, dims=1))

#     # Total N at risk
#     N = vec(sum(atrisk, dims=1))

#     # Observed Rate
#     Or = O ./ N

#     # Expected values
#     E = zeros(Float64, n_groups, ntimes)
#     for i in 1:n_groups
#         E[i, :] = Or .* vec(atrisk[i, :])
#     end

#     o = vec(sum(events, dims=2))
#     e = vec(sum(E, dims=2))

#     chi2 = sum((o .- e) .^ 2 ./ e)

#     dof = (length(o) - 1) * (length(e) - 1)
#     pval = ccdf(Chisq(dof), chi2)

#     return (o, e, size(df,1), chi2, dof, pval)
# end


