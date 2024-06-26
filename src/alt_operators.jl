#-----------------------------------------------
# functions to be used as alternative operators
# for numeric comparisons that take care of
# missing values. Comparisons involving missing
# values return "false" and any other comparisons
# return values that you would expect from
# Julia operators of ==, <, <=, >, >=.
# Broadcast operators also work. These need to be
# used as function until user-defined infix
# operators are allowed in Julia.
#-----------------------------------------------

"""
    lift(a::AbstractArray)

Takes boolean array with missing values and returns
a boolean array that replaces missing values to `false`.
"""
function lift(a::AbstractArray)
    return coalesce.(a, false)
end
