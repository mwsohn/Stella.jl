# adapted from _roc function in ROC.jl to speed up the computation

function rocdata(scores, labels)

    # thresholds are sorted list of unique values in scores
    thresholds = sort(push!(vcat(unique(scores), 1.0)), rev=true)

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
        predicted_positive = labels[mask]
        predicted_negative = labels[.!mask]
        TP[i] = sum(predicted_positive)
        TN[i] = sum(.!predicted_negative)
        FP[i] = length(predicted_positive) - TP[i]
        FN[i] = length(predicted_negative) - TN[i]
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