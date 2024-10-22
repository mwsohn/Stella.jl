function rocinput(logitmodel)
    res = response(logitmodel)
    pred = predict(logitmodel)
    mask = res .== 1.0
    return (pred[mask], pred[.!mask])
end

import Plots.plot
function plot(rocdata::Roc)
    plt = plot(rocdata.pfa, 1.0 .- rocdata.pmiss, linetype=:steppre, legend=false, xlabel="1 - specificity", ylabel="sensitivity")
    plot!(plt, collect(0:1:1), collect(0:1:1), legend=false, lc = :black)
    return plt
end

function rocplot(glmout)
    return plot(ROCAnalysis.roc(Stella.rocinput(glmout)...))
end

function roccomp(glmout1, glmout2)

    (tar1,nontar1) = rocinput(glmout1)
    (tar2,nontar2) = rocinput(glmout2)

    return ROCAnalysis.delong_test(tar1, nontar1, tar2, nontar2)
end


function auc(glmout; rocplot=false)
    rr = ROCAnalysis.roc(rocinput(glmout)...)
    println("Number of observations = ", nobs(glmout))
    println("Area under ROC curve   = ", ROCAnalysis.AUC(rr))
    if rocplot
        plt = plot(rr)
        return plt
    end
end
