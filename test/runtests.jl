using Base.Test
using Stella

# equal to
@test eq(missing, 1) == false
@test eq(1, missing) == false
@test eq(missing, missing) == false
@test eq(1, 1.0) == true
@test eq(1.0, 1.0) == true
@test eq.([1, 3, 4], [missing, 3, 1.0]) == [false, true, false]

# less than
@test lt(missing, 1) == false
@test lt(1, missing) == false
@test lt(missing, missing) == false
@test lt(1, 2.0) == true
@test lt(1.0, 2.0) == true
@test lt.([1, 3, 4], [missing, 4, 1.0]) == [false, true, false]

# less than or equal to
@test le(missing, 1) == false
@test le(1, missing) == false
@test le(missing, missing) == false
@test le(2, 1.0) == true
@test le(1, 1.0) == true
@test le(1.0, 1.0) == true
@test le.([1, 3, 4], [missing, 3, 1.0]) == [false, true, false]

# greater than
@test gt(missing, 1) == false
@test gt(1, missing) == false
@test gt(missing, missing) == false
@test gt(2.0, 1) == true
@test gt(2.0, 1.0) == true
@test gt.([1, 3, 4], [missing, 1, 1.0]) == [false, true, true]

# greater than or equal to
@test ge(missing, 1) == false
@test ge(1, missing) == false
@test ge(missing, missing) == false
@test ge(1.0, 2) == false
@test ge(2, 2.0) == true
@test ge(1.0, 1.0) == true
@test ge.([1, 3, 4], [missing, 3, 1.0]) == [false, true, true]

