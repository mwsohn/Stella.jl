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

function strpr(n::Real,d::Int)
    if d == 2
        return @sprintf("%.1f",n)
    elseif d == 3
        return @sprintf("%.2f",n)
    elseif d == 4
        return @sprintf("%.3f",n)
    elseif d == 5
        return @sprintf("%.4f",n)
    elseif d == 6
        return @sprintf("%.5f",n)
    elseif d == 7
        return @sprintf("%.6f",n)
    end
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

Base.show(io::IO,::MIME"text/plain",t::TTReturn) = show(io,t)
function show(io::IO,t::TTReturn)

    # compute widths for each column
    width,digits = _getwidths(t.colnms,t.array)

    tlen = sum(width) + 3*length(width) + 1 # 2 spaces before and after `|` + 1 boundary char

    # title
    println(io,t.title)

    # header row
    print(io,"| ",lpad(t.colnms[1],width[1]))
    for i=2:length(t.colnms)
        print(io," | ",lpad(t.colnms[i],width[i]))
    end
    print(io," |\n")

    # line
    hline(io,width)

    # Array values
    alen = t.title == "One-sample t test" ? 1 : 2
    for j=1:alen
        print(io,"| ",lpad(t.array[1][j],width[1]))
        print(io," | ",lpad(t.array[2][j],width[2]))
        for i=3:length(t.array)
            print(io," | ",lpad(strpr(t.array[i][j],digits), width[i]))
        end
        print(io," |\n")
    end

    # line
    if t.title != "One-sample t test"
        hline(io,width)
    end

    # combined - do not print if paired t test
    if t.paired == false && t.title != "One-sample t test"
        print(io,"| ",lpad("combined",width[1]))
        print(io," | ",lpad(t.array[2][3],width[2]))
        for i=3:length(t.array)
            print(io," | ",lpad(strpr(t.array[i][3],digits),width[i]))
        end
        print(io," |\n")

        # line
        hline(io,width)
    end

    # diff
    if t.title != "One-sample t test"
        print(io,"| ",lpad("diff",width[1]))
        for i=2:length(t.array)
            if t.paired == false && i in [2,4]
                print(io," | ",repeat(" ",width[i]))
            else
                print(io," | ",lpad(strpr(t.array[i][4],digits),width[i]))
            end
        end
        print(io," |\n\n")
    end

    # print ttest information
    if t.title == "One-sample t test"
        println(io,"diff = mean(",t.array[1][1],")")
        println(io,"H₀: diff = ",t.μ0)
    else
        println(io,"diff = mean(",t.array[1][1],") - mean(",t.array[1][2],")")
        println(io,"H₀: diff = 0")
    end
    if t.welch == false
        println(io,"t = ",strpr(t.t,5)," (df = ",round(Int,t.dof),")")
    else
        println(io,"t = ",strpr(t.t,5)," (df = ",strpr(t.dof,5),")")
    end
    print(io,"\n")

    tlen3rd = floor(Int,tlen/3)
    addsp = tlen % 3

    d = t.title == "One-sample t test" ? string(t.μ0) : "0"
    println(io,rpad(string("Hₐ: diff < ",d),tlen3rd),
        cpad(string("Hₐ: diff != ",d),tlen3rd),
        lpad(string("Hₐ: diff > ",d),tlen3rd+addsp)
    )
    print(io,rpad(string("Pr(T < t) = ",strpr(t.p_left,5)),tlen3rd),
        cpad(string("Pr(|T| > |t|) = ",strpr(t.p_both,5)),tlen3rd),
        lpad(string("Pr(T > t) = ",strpr(t.p_right,5)),tlen3rd+addsp)
    )
    print(io,"\n\n")
end

Base.show(io::IO,::MIME"text/plain",x::XsqReturn) = show(io,x)
function show(io::IO,x::XsqReturn)
    println(io,"Pearson χ-square (",x.dof,") = ",@sprintf("%.5f",x.chisq),"   P-value = ",@sprintf("%.5f",x.p))
end

#
# Base.show(io::IO, ::MIME"text/plain", tr::Stella.tab_return) = show(io,tr)
#
# function show(io::IO,tr::Stella.tab_return)
#
#     if length(tr.dimnames) == 1
#         show_oneway(io,NamedArray(tr.array,tr.dicts,tr.dimnames), total = total, precision = precision)
#         exit()
#     elseif length(tr.dimnames) > 2
#         error("Only up to two dimensional arrays are currently supported")
#     end
#
#     dimnames = string(tr.dimnames[1]) * " \\ " * string(tr.dimnames[2])
#     print(io,dimnames,"\n")
#
#     # temporary
#     row = col = cell = false
#
#     # total is true when row, col, or cell is true
#     if row || col || cell
#       total = true
#     end
#
#     # row names
#     rownames = collect(keys(tr.dicts[1]))
#
#     maxrowname = 5
#     for i = 1:length(rownames)
#         maxrowname = max(maxrowname,length(rownames[i]))
#     end
#
#     # column names
#     colnames = collect(keys(tr.dicts[20]))
#
#     maxcolname = 3
#     for i = 1:length(colnames)
#         maxcolname = max(maxcolname,length(colnames[i]))
#     end
#
#     # width of data columns - the same as the length of tht grand total
#     tot = sum(tr.array) # grand total
#     colwidth = length(digits(Int(floor(tot))))
#
#     # number of columns
#     ncols = length(colnames)
#
#     # number of rows
#     nrows = length(rownames)
#
#     # floating point numbers with three digits after decimal point
#     if eltype(tr.array) <: AbstractFloat
#       colwidth += 3
#     end
#
#     # determine column widths
#     colwidth = max(maxcolname,colwidth)
#
#     #---------------------------------------------------
#     # print header
#     print(io,repeat(" ",maxrowname)," |")
#
#     for i = 1:length(colnames)
#         print(io," ",lpad(string(colnames[i]),colwidth))
#     end
#
#     if total
#       print(io," | ",lpad("Total",colwidth))
#     end
#     print(io,"\n")
#
#     print(io,repeat("-",maxrowname),"-+-",repeat("-",(colwidth+1)*(length(colnames))-1))
#
#     if total
#       print(io,"-+-",repeat("-",colwidth))
#     end
#     print(io,"\n")
#     #---------------------------------------------------
#
#     # column totals
#     colsum = sum(tr.array,1)
#
#     # row totals
#     rowsum = sum(tr.array,2)
#
#     #----------------------------------------------------
#     # print values
#     for i = 1:nrows
#
#         # row name
#         print(io,rpad(string(rownames[i]),maxrowname)," |")
#
#         for j = 1:ncols
#             val = strval(tr.array[i,j])
#             print(io," ",lpad(val,colwidth))
#         end
#
#         # row total
#         if total
#           print(io," |")
#
#           val = strval(rowsum[i])
#           print(io," ",lpad(val,colwidth))
#         end
#         print(io,"\n")
#
#         # row percentages
#         if row
#           print(io,repeat(" ",maxrowname)," |")
#           for j = 1:ncols
#               val = strval(100 * tr.array[i,j] / rowsum[i],precision)
#               print(io," ",lpad(val,colwidth," "))
#           end
#
#           # row percentage
#           print(io," |")
#
#           val = strval(100.0,precision)
#           print(io," ",lpad(val,colwidth),"\n")
#         end
#
#         # column percentages
#         if col
#           print(io,repeat(" ",maxrowname)," |")
#           for j = 1:ncols
#               val = strval(100 * tr.array[i,j] / colsum[j],precision)
#               print(io," ",lpad(val,colwidth))
#           end
#
#           # column percent
#           print(" |")
#
#           val = strval(100 * rowsum[i] / tot,precision)
#           print(io," ",lpad(val,colwidth),"\n")
#         end
#
#         # column percentages
#         if cell
#           print(io,repeat(" ",maxrowname)," |")
#           for j = 1:ncols
#               val = strval(100 * tr.array[i,j] / tot,precision)
#               print(io," ",lpad(val,colwidth))
#           end
#
#           # column percent
#           print(io," |")
#
#           val = strval(100 * rowsum[i] / tot,precision)
#           print(io," ",lpad(val,colwidth),"\n")
#         end
#
#         if row || col || cell
#           print(io,repeat("-",maxrowname),"-+=",repeat("-",(colwidth+1)*(length(colnames))),"+-",repeat("-",colwidth),"\n")
#         end
#
#     end
#
#     #----------------------------------------------------
#     # Total
#     if total
#       if !(row || col || cell)
#         print(io,repeat("-",maxrowname+1),"+",repeat("-",(colwidth+1)*(length(colnames))),"-+-",repeat("-",colwidth),"\n")
#       end
#       print(rpad("Total",maxrowname," ")," |")
#
#       for i = 1:length(colnames)
#           val=strval(colsum[i])
#           print(io," ",lpad(val,colwidth," "))
#       end
#
#       # Grand total
#       val = strval(tot)
#       print(io," | ",lpad(val,colwidth," "),"\n")
#     end
#     #----------------------------------------------------
#
#     # row percentages
#     if row
#       print(io,repeat(" ",maxrowname)," |")
#       for j = 1:ncols
#           val = strval(100 * colsum[j] / tot,precision)
#           print(io," ",lpad(val,colwidth," "))
#       end
#
#       # column percent
#       val = strval(100.0,precision)
#       print(io," | ",lpad(val,colwidth," "),"\n")
#     end
#
#     # column percentages
#     if col
#       print(io,repeat(" ",maxrowname)," |")
#       for j = 1:ncols
#           val = strval(100.,precision)
#           print(io," ",lpad(val,colwidth," "))
#       end
#
#       # column percent
#       val = strval(100.0,precision)
#       print(io," | ",lpad(val,colwidth," "),"\n")
#
#     end
#
#     # column percentages
#     if cell
#       print(io,repeat(" ",maxrowname)," |")
#       for j = 1:ncols
#           val = strval(100 * colsum[j] / tot,precision)
#           print(io," ",lpad(val,colwidth," "))
#       end
#
#       # column percent
#       val = strval(100.0,precision)
#       print(io," | ",lpad(val,colwidth," "),"\n")
#     end
#
#     # p-value
#     println(io,"\nPearson χ² (",tr.dof,") = ",@sprintf("%.3f",tr.chisq)," Pr = ",@sprintf("%.3f",tr.p))
# end
#
# function show_oneway(io::IO,na::NamedArrays.NamedArray; total = false, precision::Int8 = 2)
#
#     rownames = names(na,1)
#
#     maxrowname = max(5,maximum(length.(rownames)))
#
#     # width of data columns - the same as the length of tht grand total
#     tot = sum(na) # grand total
#     colwidth = max(7,length(digits(Int(floor(tot)))))
#
#     # print header
#     print(io,rpad("Value",maxrowname," "))
#     print(io," | ",lpad("Count",colwidth," ")," ",lpad("Percent",colwidth," "),
#         " ",lpad("Cum_pct",colwidth," "),"\n")
#     print(io,repeat("-",maxrowname),"-+-",repeat("-",colwidth),"-",repeat("-",colwidth),"-",repeat("-",colwidth),"\n")
#
#     # values
#     cumtot = cumsum(na.array)
#     for i = 1:size(na.array,1)
#         str = strval(na.array[i])
#         pct = strval(100*na.array[i] / tot, precision)
#         cumpct = strval(100*cumtot[i] / tot, precision)
#         print(io,rpad(string(rownames[i]),maxrowname," "),
#             " | ",lpad(str,colwidth," "),
#             " ",lpad(pct,colwidth," "),
#             " ",lpad(cumpct,colwidth," "),"\n")
#     end
#
#     # total
#     if total
#         print(io,repeat("-",maxrowname),"-+-",repeat("-",colwidth),"-",repeat("-",colwidth),"-",repeat("-",colwidth),"\n")
#         str = strval(tot)
#         pct = strval(100.00, precision)
#         cumpct = strval(100.00, precision)
#         print(io,rpad("Total",maxrowname," "),
#             " | ",lpad(str,colwidth," "),
#             " ",lpad(pct,colwidth," "),
#             " ",lpad(cumpct,colwidth," "),"\n")
#     end
# end
#
# Base.show(io::IO,::MIME"text/plain",pr::pwcorr_return) = show(io,pr)
#
# function show(io::IO,pr::pwcorr_return; width::Int8 = 9, precision::Int8 = 3, p = false, N = false)
#
#     ncol = size(pr.N,1)
#
#     for i = 1:ncol
#         print(io,prepend_spaces(pr.colnames[i],width))
#         if i < ncol
#             print(io," ")
#         end
#     end
#     print(io,"\n")
#
#     for i = 1:ncol
#         print(io,prepend_spaces(pr.colnames[i],width)," ")
#         for j=1:i
#             printf(io,"%9.4f ",bc[i,j])
#             if i == j
#                 print(io,"\n")
#             end
#         end
#         if p == true
#             print(io,repeat(" ",width+1))
#             for j=1:i
#                 printf(io,"%9.4f ",bp[i,j])
#                 if i == j
#                     print(io,"\n")
#                 end
#             end
#         end
#     end
# end
# #
# # function show(io::IO,tr::Stella.tab_return; row=false, col=false, cell=false, total=false, precision::Int8 = 2)
# #
# #     if length(tr.dimnames) == 1
# #         show_oneway(io,NamedArray(tr.array,tr.dicts,tr.dimnames), total = total, precision = precision)
# #         exit()
# #     elseif length(tr.dimnames) > 2
# #         error("Only up to two dimensional arrays are currently supported")
# #     end
# #
# #     dimnames = string(tr.dimnames[1]) * " \\ " * string(tr.dimnames[2])
# #     print(io,dimnames,"\n")
# #
# #     # total is true when row, col, or cell is true
# #     if row || col || cell
# #       total = true
# #     end
# #
# #     # row names
# #     rownames = collect(keys(tr.dicts[1]))
# #
# #     maxrowname = 5
# #     for i = 1:length(rownames)
# #         maxrowname = max(maxrowname,length(rownames[i]))
# #     end
# #
# #     # column names
# #     colnames = collect(keys(tr.dicts[20]))
# #
# #     maxcolname = 3
# #     for i = 1:length(colnames)
# #         maxcolname = max(maxcolname,length(colnames[i]))
# #     end
# #
# #     # width of data columns - the same as the length of tht grand total
# #     tot = sum(tr.array) # grand total
# #     colwidth = length(digits(Int(floor(tot))))
# #
# #     # number of columns
# #     ncols = length(colnames)
# #
# #     # number of rows
# #     nrows = length(rownames)
# #
# #     # floating point numbers with three digits after decimal point
# #     if eltype(tr.array) <: AbstractFloat
# #       colwidth += 3
# #     end
# #
# #     # determine column widths
# #     colwidth = max(maxcolname,colwidth)
# #
# #     #---------------------------------------------------
# #     # print header
# #     print(io,repeat(" ",maxrowname)," |")
# #
# #     for i = 1:length(colnames)
# #         print(io," ",lpad(string(colnames[i]),colwidth))
# #     end
# #
# #     if total
# #       print(io," | ",lpad("Total",colwidth))
# #     end
# #     print(io,"\n")
# #
# #     print(io,repeat("-",maxrowname),"-+-",repeat("-",(colwidth+1)*(length(colnames))-1))
# #
# #     if total
# #       print(io,"-+-",repeat("-",colwidth))
# #     end
# #     print(io,"\n")
# #     #---------------------------------------------------
# #
# #     # column totals
# #     colsum = sum(tr.array,1)
# #
# #     # row totals
# #     rowsum = sum(tr.array,2)
# #
# #     #----------------------------------------------------
# #     # print values
# #     for i = 1:nrows
# #
# #         # row name
# #         print(io,rpad(string(rownames[i]),maxrowname)," |")
# #
# #         for j = 1:ncols
# #             val = strval(tr.array[i,j])
# #             print(io," ",lpad(val,colwidth))
# #         end
# #
# #         # row total
# #         if total
# #           print(io," |")
# #
# #           val = strval(rowsum[i])
# #           print(io," ",lpad(val,colwidth))
# #         end
# #         print(io,"\n")
# #
# #         # row percentages
# #         if row
# #           print(io,repeat(" ",maxrowname)," |")
# #           for j = 1:ncols
# #               val = strval(100 * tr.array[i,j] / rowsum[i],precision)
# #               print(io," ",lpad(val,colwidth," "))
# #           end
# #
# #           # row percentage
# #           print(io," |")
# #
# #           val = strval(100.0,precision)
# #           print(io," ",lpad(val,colwidth),"\n")
# #         end
# #
# #         # column percentages
# #         if col
# #           print(io,repeat(" ",maxrowname)," |")
# #           for j = 1:ncols
# #               val = strval(100 * tr.array[i,j] / colsum[j],precision)
# #               print(io," ",lpad(val,colwidth))
# #           end
# #
# #           # column percent
# #           print(" |")
# #
# #           val = strval(100 * rowsum[i] / tot,precision)
# #           print(io," ",lpad(val,colwidth),"\n")
# #         end
# #
# #         # column percentages
# #         if cell
# #           print(io,repeat(" ",maxrowname)," |")
# #           for j = 1:ncols
# #               val = strval(100 * tr.array[i,j] / tot,precision)
# #               print(io," ",lpad(val,colwidth))
# #           end
# #
# #           # column percent
# #           print(io," |")
# #
# #           val = strval(100 * rowsum[i] / tot,precision)
# #           print(io," ",lpad(val,colwidth),"\n")
# #         end
# #
# #         if row || col || cell
# #           print(io,repeat("-",maxrowname),"-+=",repeat("-",(colwidth+1)*(length(colnames))),"+-",repeat("-",colwidth),"\n")
# #         end
# #
# #     end
# #
# #     #----------------------------------------------------
# #     # Total
# #     if total
# #       if !(row || col || cell)
# #         print(io,repeat("-",maxrowname+1),"+",repeat("-",(colwidth+1)*(length(colnames))),"-+-",repeat("-",colwidth),"\n")
# #       end
# #       print(rpad("Total",maxrowname," ")," |")
# #
# #       for i = 1:length(colnames)
# #           val=strval(colsum[i])
# #           print(io," ",lpad(val,colwidth," "))
# #       end
# #
# #       # Grand total
# #       val = strval(tot)
# #       print(io," | ",lpad(val,colwidth," "),"\n")
# #     end
# #     #----------------------------------------------------
# #
# #     # row percentages
# #     if row
# #       print(io,repeat(" ",maxrowname)," |")
# #       for j = 1:ncols
# #           val = strval(100 * colsum[j] / tot,precision)
# #           print(io," ",lpad(val,colwidth," "))
# #       end
# #
# #       # column percent
# #       val = strval(100.0,precision)
# #       print(io," | ",lpad(val,colwidth," "),"\n")
# #     end
# #
# #     # column percentages
# #     if col
# #       print(io,repeat(" ",maxrowname)," |")
# #       for j = 1:ncols
# #           val = strval(100.,precision)
# #           print(io," ",lpad(val,colwidth," "))
# #       end
# #
# #       # column percent
# #       val = strval(100.0,precision)
# #       print(io," | ",lpad(val,colwidth," "),"\n")
# #
# #     end
# #
# #     # column percentages
# #     if cell
# #       print(io,repeat(" ",maxrowname)," |")
# #       for j = 1:ncols
# #           val = strval(100 * colsum[j] / tot,precision)
# #           print(io," ",lpad(val,colwidth," "))
# #       end
# #
# #       # column percent
# #       val = strval(100.0,precision)
# #       print(io," | ",lpad(val,colwidth," "),"\n")
# #     end
# #
# #     # p-value
# #     println(io,"\nPearson χ² (",tr.dof,") = ",@sprintf("%.3f",tr.chisq)," Pr = ",@sprintf("%.3f",tr.p))
# # end
# #
# # function show_oneway(io::IO,na::NamedArrays.NamedArray; total = false, precision::Int8 = 2)
# #
# #     rownames = names(na,1)
# #
# #     maxrowname = max(5,maximum(length.(rownames)))
# #
# #     # width of data columns - the same as the length of tht grand total
# #     tot = sum(na) # grand total
# #     colwidth = max(7,length(digits(Int(floor(tot)))))
# #
# #     # print header
# #     print(io,rpad("Value",maxrowname," "))
# #     print(io," | ",lpad("Count",colwidth," ")," ",lpad("Percent",colwidth," "),
# #         " ",lpad("Cum_pct",colwidth," "),"\n")
# #     print(io,repeat("-",maxrowname),"-+-",repeat("-",colwidth),"-",repeat("-",colwidth),"-",repeat("-",colwidth),"\n")
# #
# #     # values
# #     cumtot = cumsum(na.array)
# #     for i = 1:size(na.array,1)
# #         str = strval(na.array[i])
# #         pct = strval(100*na.array[i] / tot, precision)
# #         cumpct = strval(100*cumtot[i] / tot, precision)
# #         print(io,rpad(string(rownames[i]),maxrowname," "),
# #             " | ",lpad(str,colwidth," "),
# #             " ",lpad(pct,colwidth," "),
# #             " ",lpad(cumpct,colwidth," "),"\n")
# #     end
# #
# #     # total
# #     if total
# #         print(io,repeat("-",maxrowname),"-+-",repeat("-",colwidth),"-",repeat("-",colwidth),"-",repeat("-",colwidth),"\n")
# #         str = strval(tot)
# #         pct = strval(100.00, precision)
# #         cumpct = strval(100.00, precision)
# #         print(io,rpad("Total",maxrowname," "),
# #             " | ",lpad(str,colwidth," "),
# #             " ",lpad(pct,colwidth," "),
# #             " ",lpad(cumpct,colwidth," "),"\n")
# #     end
# # end
# #
# # Base.show(io::IO,::MIME"text/plain",pr::pwcorr_return) = show(io,pr)
# #
# # function show(io::IO,pr::pwcorr_return; width::Int8 = 9, precision::Int8 = 3, p = false, N = false)
# #
# #     ncol = size(pr.N,1)
# #
# #     for i = 1:ncol
# #         print(io,prepend_spaces(pr.colnames[i],width))
# #         if i < ncol
# #             print(io," ")
# #         end
# #     end
# #     print(io,"\n")
# #
# #     for i = 1:ncol
# #         print(io,prepend_spaces(pr.colnames[i],width)," ")
# #         for j=1:i
# #             printf(io,"%9.4f ",bc[i,j])
# #             if i == j
# #                 print(io,"\n")
# #             end
# #         end
# #         if p == true
# #             print(io,repeat(" ",width+1))
# #             for j=1:i
# #                 printf(io,"%9.4f ",bp[i,j])
# #                 if i == j
# #                     print(io,"\n")
# #                 end
# #             end
# #         end
# #     end
# # end
