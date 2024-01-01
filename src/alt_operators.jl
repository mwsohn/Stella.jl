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
    lift(a::AbstractArray{Union{Bool,Missing}})

Takes boolean array with missing values and returns
a boolean array that replaces missing values to `false`.
"""
function lift(a::AbstractArray{Union{Bool,Missing}})
    return coalesce.(a, false)
end

function eq(a::Missing, b::Real)
    return false
end
function eq(a::Real, b::Missing)
    return false
end
function eq(a::Missing, b::Missing)
    return false
end
function eq(a::Real, b::Real)
    return a == b
end
function eq(a::Missing, b::AbstractString)
    return false
end
function eq(a::AbstractString, b::Missing)
    return false
end
function eq(a::AbstractString, b::AbstractString)
    return a == b
end

function lt(a::Missing, b::Real)
    return false
end
function lt(a::Real, b::Missing)
    return false
end
function lt(a::Missing, b::Missing)
    return false
end
function lt(a::Real, b::Real)
    return a < b
end

function le(a::Missing, b::Real)
    return false
end
function le(a::Real, b::Missing)
    return false
end
function le(a::Missing, b::Missing)
    return false
end
function le(a::Real, b::Real)
    return a <= b
end

function gt(a::Missing, b::Real)
    return false
end
function gt(a::Real, b::Missing)
    return false
end
function gt(a::Missing, b::Missing)
    return false
end
function gt(a::Real, b::Real)
    return a > b
end

function ge(a::Missing, b::Real)
    return false
end
function ge(a::Real, b::Missing)
    return false
end
function ge(a::Missing, b::Missing)
    return false
end
function ge(a::Real, b::Real)
    return a >= b
end
