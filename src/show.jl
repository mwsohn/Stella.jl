import Base.show

function hline(io,width::Vector)
    print(io,"├")
    for i=1:length(width)
        print(io,repeat("─",width[i]+2))
        if i < length(width)
            print(io,"┼")
        end
    end
    print(io,"┤\n")
end

function cpad(str::String,w::Int)
    len = length(str)
    if len > w
        return str[1:w]
    end

    d = w - len
    lsp = rsp = floor(Int,d/2)
    if d % 2 != 0
        rsp += 1
    end
    return string(repeat(" ",lsp),str,repeat(" ",rsp))
end

function _getwidths(hdr,varray)
    width = length.(hdr)
    width[1] = max(width[1],maximum(length.(string.(varray[1]))),8) # 8 is for "combined"
    width[2] = max(width[2],maximum(length(string.(varray[2]))),6) # 6 is minimum width for any column
    mlen = max(length(string(maximum(maximum.(round.(Int,varray[3:end]))))),
        length(string(minimum(minimum.(round.(Int,varray[3:end]))))))
    digits = mlen > 5 ? 3 : mlen > 4 ? 4 : mlen > 3 ? 4 : 5 # minimum two decimal places
    for i=3:length(varray) # varray needs to be vector of arrays
        mlen = maximum(length.(string.(round.(Int,varray[i]))))
        width[i] = max(width[i],8,mlen + digits)
    end
    return (width,digits)
end

# Base.show(io::IO,::MIME"text/plain",t::TTReturn) = show(io,t)
# function Base.show(io::IO,t::TTReturn)

#     # compute widths for each column
#     width,digits = _getwidths(t.colnms,t.array)

#     tlen = sum(width) + 3*length(width) + 1 # 2 spaces before and after `|` + 1 boundary char

#     # title
#     println(io,t.title)

#     # header row
#     print(io,"| ",lpad(t.colnms[1],width[1]))
#     for i=2:length(t.colnms)
#         print(io," | ",lpad(t.colnms[i],width[i]))
#     end
#     print(io," |\n")

#     # line
#     hline(io,width)

#     # Array values
#     alen = t.title == "One-sample t test" ? 1 : 2
#     for j=1:alen
#         print(io,"| ",lpad(t.array[1][j],width[1]))
#         print(io," | ",lpad(t.array[2][j],width[2]))
#         for i=3:length(t.array)
#             print(io," | ",lpad(strpr(t.array[i][j],digits), width[i]))
#         end
#         print(io," |\n")
#     end

#     # line
#     if t.title != "One-sample t test"
#         hline(io,width)
#     end

#     # combined - do not print if paired t test
#     if t.paired == false && t.title != "One-sample t test"
#         print(io,"| ",lpad("combined",width[1]))
#         print(io," | ",lpad(t.array[2][3],width[2]))
#         for i=3:length(t.array)
#             print(io," | ",lpad(strpr(t.array[i][3],digits),width[i]))
#         end
#         print(io," |\n")

#         # line
#         hline(io,width)
#     end

#     # diff
#     if t.title != "One-sample t test"
#         print(io,"| ",lpad("diff",width[1]))
#         for i=2:length(t.array)
#             if t.paired == false && i in [2,4]
#                 print(io," | ",repeat(" ",width[i]))
#             else
#                 print(io," | ",lpad(strpr(t.array[i][4],digits),width[i]))
#             end
#         end
#         print(io," |\n\n")
#     end

#     # print ttest information
#     if t.title == "One-sample t test"
#         println(io,"diff = mean(",t.array[1][1],")")
#         println(io,"H₀: diff = ",t.μ0)
#     else
#         println(io,"diff = mean(",t.array[1][1],") - mean(",t.array[1][2],")")
#         println(io,"H₀: diff = 0")
#     end
#     if t.welch == false
#         println(io,"t = ",strpr(t.t,5)," (df = ",round(Int,t.dof),")")
#     else
#         println(io,"t = ",strpr(t.t,5)," (df = ",strpr(t.dof,5),")")
#     end
#     print(io,"\n")

#     tlen3rd = floor(Int,tlen/3)
#     addsp = tlen % 3

#     d = t.title == "One-sample t test" ? string(t.μ0) : "0"
#     println(io,rpad(string("Hₐ: diff < ",d),tlen3rd),
#         cpad(string("Hₐ: diff != ",d),tlen3rd),
#         lpad(string("Hₐ: diff > ",d),tlen3rd+addsp)
#     )
#     print(io,rpad(string("Pr(T < t) = ",strpr(t.p_left,5)),tlen3rd),
#         cpad(string("Pr(|T| > |t|) = ",strpr(t.p_both,5)),tlen3rd),
#         lpad(string("Pr(T > t) = ",strpr(t.p_right,5)),tlen3rd+addsp)
#     )
#     print(io,"\n\n")
# end

# Base.show(io::IO,::MIME"text/plain",x::XsqResult) = show(io,x)
function Base.show(io::IO,x::XsqResult)
    println(io,"Pearson χ-square (",x.dof,") = ",@sprintf("%.5f",x.chisq),"   P-value = ",@sprintf("%.5f",x.p))
end

