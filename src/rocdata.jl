# adapted from _roc function in ROC.jl to speed up the computation

function rocdata(scores, labels)
    thresholds = sort(push!(unique(scores), 1), rev=true)
    labels = labels .== 1.0 # turn it into a boolean vector
    P = sum(labels)        # positives
    N = length(labels) - P # negatives
    n_thresholds = length(thresholds)

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
        # TP[i] = TPi
        # TN[i] = TNi
        FP[i] = length(predicted_positive) - TP[i]
        FN[i] = length(predicted_negative) - TN[i]
        FPR[i] = FP[i] / (FP[i] + TN[i])
        TPR[i] = TP[i] / (TP[i] + FN[i])
    end
    ROCData{eltype(thresholds)}(thresholds, P, N, TP, TN, FP, FN, FPR, TPR)
end

auc(scores, labels) = ROC.AUC(rocdata(scores,labels))