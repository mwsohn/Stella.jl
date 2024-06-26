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
eq(missing, b) = false
eq(a, missing) = false
eq(missing, missing) = false
eq(a::Real, b::Real) = ==(a, b)
eq(missing, b::AbstractString) = false
eq(missing, b::Union{Missing,AbstractString}) = false
eq(a::AbstractString, missing) = false
eq(a::Union{Missing,AbstractString}, missing) = false
eq(a::AbstractString, b::AbstractString) = a == b
eq(a::Union{Missing,AbstractString}, b::Union{Missing,AbstractString}) = a == b

lt(a::Missing, b::Real) = false
lt(a::Real, b::Missing) = false
lt(a::Missing, b::Missing) = false
lt(a::Real, b::Real) = a < b

le(a::Missing, b::Real) = false
le(a::Real, b::Missing) = false
le(a::Missing, b::Missing) = false
le(a::Real, b::Real) = a <= b

gt(a::Missing, b::Real) = false
gt(a::Real, b::Missing) = false
gt(a::Missing, b::Missing) = false
gt(a::Real, b::Real) = a > b

ge(a::Missing, b::Real) = false
ge(a::Real, b::Missing) = false
ge(a::Missing, b::Missing) = false
ge(a::Real, b::Real) = a >= b

