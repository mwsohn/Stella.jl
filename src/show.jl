import Base.show

Base.show(io::IO, ::MIME"text/plain", t::Stella.tab_return) = show(io,t)

function show(io::IO,tr::Stella.tab_return; row=false, col=false, cell=false, total=false, precision::Int8 = 2)

    if length(tr.dimnames) == 1
        show_oneway(io,NamedArray(tr.array,tr.dicts,tr.dimnames), total = total, precision = precision)
        exit()
    elseif length(tr.dimnames) > 2
        error("Only up to two dimensional arrays are currently supported")
    end

    dimnames = string(tr.dimnames[1]) * " \\ " * string(tr.dimnames[2])
    print(io,dimnames,"\n")

    # total is true when row, col, or cell is true
    if row || col || cell
      total = true
    end

    # row names
    rownames = collect(keys(tr.dicts[1]))

    maxrowname = 5
    for i = 1:length(rownames)
        maxrowname = max(maxrowname,length(rownames[i]))
    end

    # column names
    colnames = collect(keys(tr.dicts[20]))

    maxcolname = 3
    for i = 1:length(colnames)
        maxcolname = max(maxcolname,length(colnames[i]))
    end

    # width of data columns - the same as the length of tht grand total
    tot = sum(tr.array) # grand total
    colwidth = length(digits(Int(floor(tot))))

    # number of columns
    ncols = length(colnames)

    # number of rows
    nrows = length(rownames)

    # floating point numbers with three digits after decimal point
    if eltype(tr.array) <: AbstractFloat
      colwidth += 3
    end

    # determine column widths
    colwidth = max(maxcolname,colwidth)

    #---------------------------------------------------
    # print header
    print(io,repeat(" ",maxrowname)," |")

    for i = 1:length(colnames)
        print(io," ",lpad(string(colnames[i]),colwidth))
    end

    if total
      print(io," | ",lpad("Total",colwidth))
    end
    print(io,"\n")

    print(io,repeat("-",maxrowname),"-+-",repeat("-",(colwidth+1)*(length(colnames))-1))

    if total
      print(io,"-+-",repeat("-",colwidth))
    end
    print(io,"\n")
    #---------------------------------------------------

    # column totals
    colsum = sum(tr.array,1)

    # row totals
    rowsum = sum(tr.array,2)

    #----------------------------------------------------
    # print values
    for i = 1:nrows

        # row name
        print(io,rpad(string(rownames[i]),maxrowname)," |")

        for j = 1:ncols
            val = strval(tr.array[i,j])
            print(io," ",lpad(val,colwidth))
        end

        # row total
        if total
          print(io," |")

          val = strval(rowsum[i])
          print(io," ",lpad(val,colwidth))
        end
        print(io,"\n")

        # row percentages
        if row
          print(io,repeat(" ",maxrowname)," |")
          for j = 1:ncols
              val = strval(100 * tr.array[i,j] / rowsum[i],precision)
              print(io," ",lpad(val,colwidth," "))
          end

          # row percentage
          print(io," |")

          val = strval(100.0,precision)
          print(io," ",lpad(val,colwidth),"\n")
        end

        # column percentages
        if col
          print(io,repeat(" ",maxrowname)," |")
          for j = 1:ncols
              val = strval(100 * tr.array[i,j] / colsum[j],precision)
              print(io," ",lpad(val,colwidth))
          end

          # column percent
          print(" |")

          val = strval(100 * rowsum[i] / tot,precision)
          print(io," ",lpad(val,colwidth),"\n")
        end

        # column percentages
        if cell
          print(io,repeat(" ",maxrowname)," |")
          for j = 1:ncols
              val = strval(100 * tr.array[i,j] / tot,precision)
              print(io," ",lpad(val,colwidth))
          end

          # column percent
          print(io," |")

          val = strval(100 * rowsum[i] / tot,precision)
          print(io," ",lpad(val,colwidth),"\n")
        end

        if row || col || cell
          print(io,repeat("-",maxrowname),"-+=",repeat("-",(colwidth+1)*(length(colnames))),"+-",repeat("-",colwidth),"\n")
        end

    end

    #----------------------------------------------------
    # Total
    if total
      if !(row || col || cell)
        print(io,repeat("-",maxrowname+1),"+",repeat("-",(colwidth+1)*(length(colnames))),"-+-",repeat("-",colwidth),"\n")
      end
      print(rpad("Total",maxrowname," ")," |")

      for i = 1:length(colnames)
          val=strval(colsum[i])
          print(io," ",lpad(val,colwidth," "))
      end

      # Grand total
      val = strval(tot)
      print(io," | ",lpad(val,colwidth," "),"\n")
    end
    #----------------------------------------------------

    # row percentages
    if row
      print(io,repeat(" ",maxrowname)," |")
      for j = 1:ncols
          val = strval(100 * colsum[j] / tot,precision)
          print(io," ",lpad(val,colwidth," "))
      end

      # column percent
      val = strval(100.0,precision)
      print(io," | ",lpad(val,colwidth," "),"\n")
    end

    # column percentages
    if col
      print(io,repeat(" ",maxrowname)," |")
      for j = 1:ncols
          val = strval(100.,precision)
          print(io," ",lpad(val,colwidth," "))
      end

      # column percent
      val = strval(100.0,precision)
      print(io," | ",lpad(val,colwidth," "),"\n")

    end

    # column percentages
    if cell
      print(io,repeat(" ",maxrowname)," |")
      for j = 1:ncols
          val = strval(100 * colsum[j] / tot,precision)
          print(io," ",lpad(val,colwidth," "))
      end

      # column percent
      val = strval(100.0,precision)
      print(io," | ",lpad(val,colwidth," "),"\n")
    end

    # p-value
    println(io,"\nPearson χ² (",tr.dof,") = ",@sprintf("%.3f",tr.chisq)," Pr = ",@sprintf("%.3f",tr.p))
end

function show_oneway(io::IO,na::NamedArrays.NamedArray; total = false, precision::Int8 = 2)

    rownames = names(na,1)

    maxrowname = max(5,maximum(length.(rownames)))

    # width of data columns - the same as the length of tht grand total
    tot = sum(na) # grand total
    colwidth = max(7,length(digits(Int(floor(tot)))))

    # print header
    print(io,rpad("Value",maxrowname," "))
    print(io," | ",lpad("Count",colwidth," ")," ",lpad("Percent",colwidth," "),
        " ",lpad("Cum_pct",colwidth," "),"\n")
    print(io,repeat("-",maxrowname),"-+-",repeat("-",colwidth),"-",repeat("-",colwidth),"-",repeat("-",colwidth),"\n")

    # values
    cumtot = cumsum(na.array)
    for i = 1:size(na.array,1)
        str = strval(na.array[i])
        pct = strval(100*na.array[i] / tot, precision)
        cumpct = strval(100*cumtot[i] / tot, precision)
        print(io,rpad(string(rownames[i]),maxrowname," "),
            " | ",lpad(str,colwidth," "),
            " ",lpad(pct,colwidth," "),
            " ",lpad(cumpct,colwidth," "),"\n")
    end

    # total
    if total
        print(io,repeat("-",maxrowname),"-+-",repeat("-",colwidth),"-",repeat("-",colwidth),"-",repeat("-",colwidth),"\n")
        str = strval(tot)
        pct = strval(100.00, precision)
        cumpct = strval(100.00, precision)
        print(io,rpad("Total",maxrowname," "),
            " | ",lpad(str,colwidth," "),
            " ",lpad(pct,colwidth," "),
            " ",lpad(cumpct,colwidth," "),"\n")
    end
end
