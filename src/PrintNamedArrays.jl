# using NamedArrays, Distributions
#
# import Base.print
#
# function strval(val::AbstractFloat)
#   return @sprintf("%.2f",val)
# end
#
# function strval(val::AbstractFloat, decimals::Int)
#     if decimals == 1
#         return @sprintf("%.1f",val)
#     elseif decimals == 2
#         return @sprintf("%.2f",val)
#     elseif decimals == 3
#         return @sprintf("%.3f",val)
#     elseif decimals == 4
#         return @sprintf("%.4f",val)
#     elseif decimals == 5
#         return @sprintf("%.5f",val)
#     elseif decimals == 6
#         return @sprintf("%.6f",val)
#     else
#         return @sprintf("%f",val)
#     end
# end
#
# function strval(val::Integer)
#   return @sprintf("%.0d",val)
# end
#
# function getdictval(dt::Dict,val)
#     return haskey(dt,val) ? dt[val] : val
# end
#
# function print(tr::tab_return; row=false, col=false, cell=false, total=false, precision::Int8 = 2)
#
#     # named array
#     na = tr.na
#
#     if ndims(na.array) == 1
#         return print_oneway(na, total = total, precision = precision)
#     elseif length(na.dimnames) > 2
#         error("Only up to two dimensional arrays are currently supported")
#     end
#
#     dimnames = string(na.dimnames[1]) * " \\ " * string(na.dimnames[2])
#     print(dimnames,"\n")
#
#     # total is true when row, col, or cell is true
#     if row || col || cell
#       total = true
#     end
#
#     # row names
#     rownames = names(na,1)
#
#     maxrowname = 5
#     for i = 1:length(rownames)
#         maxrowname = max(maxrowname,length(rownames[i]))
#     end
#
#     # column names
#     colnames = names(na,2)
#
#     maxcolname = 3
#     for i = 1:length(colnames)
#         maxcolname = max(maxcolname,length(colnames[i]))
#     end
#
#     # width of data columns - the same as the length of tht grand total
#     tot = sum(na) # grand total
#     colwidth = length(digits(Int(floor(tot))))
#
#     # number of columns
#     ncols = length(colnames)
#
#     # number of rows
#     nrows = length(rownames)
#
#     # floating point numbers with three digits after decimal point
#     if eltype(na.array) <: AbstractFloat
#       colwidth += 3
#     end
#
#     # determine column widths
#     colwidth = max(maxcolname,colwidth)
#
#     #---------------------------------------------------
#     # print header
#     print(repeat(" ",maxrowname)," |")
#
#     for i = 1:length(colnames)
#         print(" ",lpad(string(colnames[i]),colwidth," "))
#     end
#
#     if total
#       print(" | ",lpad("Total",colwidth," "))
#     end
#     print("\n")
#
#     print(repeat("-",maxrowname),"-+-",repeat("-",(colwidth+1)*(length(colnames))-1))
#
#     if total
#       print("-+-",repeat("-",colwidth))
#     end
#     print("\n")
#     #---------------------------------------------------
#
#     # column totals
#     colsum = sum(na.array,1)
#
#     # row totals
#     rowsum = sum(na.array,2)
#
#     #----------------------------------------------------
#     # print values
#     for i = 1:nrows
#
#         # row name
#         print(rpad(string(rownames[i]),maxrowname," ")," |")
#
#         for j = 1:ncols
#             val = strval(na.array[i,j])
#             print(" ",lpad(val,colwidth," "))
#         end
#
#         # row total
#         if total
#           print(" |")
#
#           val = strval(rowsum[i])
#           print(" ",lpad(val,colwidth," "))
#         end
#         print("\n")
#
#         # row percentages
#         if row
#           print(repeat(" ",maxrowname)," |")
#           for j = 1:ncols
#               val = strval(100 * na.array[i,j] / rowsum[i],precision)
#               print(" ",lpad(val,colwidth," "))
#           end
#
#           # row percentage
#           print(" |")
#
#           val = strval(100.0,precision)
#           print(" ",lpad(val,colwidth," "),"\n")
#         end
#
#         # column percentages
#         if col
#           print(repeat(" ",maxrowname)," |")
#           for j = 1:ncols
#               val = strval(100 * na.array[i,j] / colsum[j],precision)
#               print(" ",lpad(val,colwidth," "))
#           end
#
#           # column percent
#           print(" |")
#
#           val = strval(100 * rowsum[i] / tot,precision)
#           print(" ",lpad(val,colwidth," "),"\n")
#         end
#
#         # column percentages
#         if cell
#           print(repeat(" ",maxrowname)," |")
#           for j = 1:ncols
#               val = strval(100 * na.array[i,j] / tot,precision)
#               print(" ",lpad(val,colwidth," "))
#           end
#
#           # column percent
#           print(" |")
#
#           val = strval(100 * rowsum[i] / tot,precision)
#           print(" ",lpad(val,colwidth," "),"\n")
#         end
#
#         if row || col || cell
#           print(repeat("-",maxrowname),"-+=",repeat("-",(colwidth+1)*(length(colnames))),"+-",repeat("-",colwidth),"\n")
#         end
#
#     end
#
#     #----------------------------------------------------
#     # Total
#     if total
#       if !(row || col || cell)
#         print(repeat("-",maxrowname+1),"+",repeat("-",(colwidth+1)*(length(colnames))),"-+-",repeat("-",colwidth),"\n")
#       end
#       print(rpad("Total",maxrowname," ")," |")
#
#       for i = 1:length(colnames)
#           val=strval(colsum[i])
#           print(" ",lpad(val,colwidth," "))
#       end
#
#       # Grand total
#       val = strval(tot)
#       print(" | ",lpad(val,colwidth," "),"\n")
#     end
#     #----------------------------------------------------
#
#     # row percentages
#     if row
#       print(repeat(" ",maxrowname)," |")
#       for j = 1:ncols
#           val = strval(100 * colsum[j] / tot,precision)
#           print(" ",lpad(val,colwidth," "))
#       end
#
#       # column percent
#       val = strval(100.0,precision)
#       print(" | ",lpad(val,colwidth," "),"\n")
#     end
#
#     # column percentages
#     if col
#       print(repeat(" ",maxrowname)," |")
#       for j = 1:ncols
#           val = strval(100.,precision)
#           print(" ",lpad(val,colwidth," "))
#       end
#
#       # column percent
#       val = strval(100.0,precision)
#       print(" | ",lpad(val,colwidth," "),"\n")
#
#     end
#
#     # column percentages
#     if cell
#       print(repeat(" ",maxrowname)," |")
#       for j = 1:ncols
#           val = strval(100 * colsum[j] / tot,precision)
#           print(" ",lpad(val,colwidth," "))
#       end
#
#       # column percent
#       val = strval(100.0,precision)
#       print(" | ",lpad(val,colwidth," "),"\n")
#     end
#
#     # p-value
#     println("\nPearson χ² (",tr.dof,") = ",@sprintf("%.3f",tr.chisq)," Pr = ",@sprintf("%.3f",tr.p))
# end
#
# function print_oneway(na::NamedArray; total = false, precision::Int8 = 2)
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
#     print(rpad("Value",maxrowname," "))
#     print(" | ",lpad("Count",colwidth," ")," ",lpad("Percent",colwidth," "),
#         " ",lpad("Cum_pct",colwidth," "),"\n")
#     print(repeat("-",maxrowname),"-+-",repeat("-",colwidth),"-",repeat("-",colwidth),"-",repeat("-",colwidth),"\n")
#
#     # values
#     cumtot = cumsum(na.array)
#     for i = 1:size(na.array,1)
#         str = strval(na.array[i])
#         pct = strval(100*na.array[i] / tot, precision)
#         cumpct = strval(100*cumtot[i] / tot, precision)
#         print(rpad(string(rownames[i]),maxrowname," "),
#             " | ",lpad(str,colwidth," "),
#             " ",lpad(pct,colwidth," "),
#             " ",lpad(cumpct,colwidth," "),"\n")
#     end
#
#     # total
#     if total
#         print(repeat("-",maxrowname),"-+-",repeat("-",colwidth),"-",repeat("-",colwidth),"-",repeat("-",colwidth),"\n")
#         str = strval(tot)
#         pct = strval(100.00, precision)
#         cumpct = strval(100.00, precision)
#         print(rpad("Total",maxrowname," "),
#             " | ",lpad(str,colwidth," "),
#             " ",lpad(pct,colwidth," "),
#             " ",lpad(cumpct,colwidth," "),"\n")
#     end
# end
