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

# alternative operators
==(a::Missing, b::Real) = false
==(a::Real, b::Missing) = false
==(a::Missing, b::Missing) = false
==(a::Real, b::Real) = a == b
==(a::Missing, b::AbstractString) = false
==(a::AbstractString, b::Missing) = false
==(a::AbstractString, b::AbstractString) = a == b

<(a::Missing, b::Real) = false
<(a::Real, b::Missing) = false
<(a::Missing, b::Missing) = false
<(a::Real, b::Real) = a < b

<=(a::Missing, b::Real) = false
<=(a::Real, b::Missing) = false
<=(a::Missing, b::Missing) = false
<=(a::Real, b::Real) = a <= b

>(a::Missing, b::Real) = false
>(a::Real, b::Missing) = false
>(a::Missing, b::Missing) = false
>(a::Real, b::Real) = a > b

>=(a::Missing, b::Real) = false
>=(a::Real, b::Missing) = false
>=(a::Missing, b::Missing) = false
>=(a::Real, b::Real) = a >= b

