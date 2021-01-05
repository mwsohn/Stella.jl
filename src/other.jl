"""
    identify_condition(df::AbstractDataFrame,cols::Vector{Symbol},codes::Vector{String})

produces a vector of Boolean `true` or `false` depending on whether columns `cols` contain
any of the string values in `codes`. This function is designed to easily identify comorbid
conditions in inpatient or outpatient records that are coded in ICD-9 or ICD-10 codes.
The values in `cols` are compared to the values (`code`) in `codes` only up to the length of
the `code`. For example, when detecting "diabetes", one may use "250.xx" in ICD-9,
which can be specified as identify_condition(df,[:diag1,:diag2,:diag3,:diag4,:diag5],["250"]).
All values that start with "250" in :diag1 - :diag4 (e.g., "25010") will match.
"""
function identify_condition(ip::AbstractDataFrame,vicd::Vector{Symbol},codes::Vector{String})

    retvec = zeros(Bool,nrow(ip))

    for i in 1:nrow(ip)
        for v in vicd
            if ismissing(ip[i,v])
                continue
            end

            vlen=length(ip[i,v])
            for c in codes
                clen = length(c)
                if vlen >= clen && c == ip[i,v][1:length(c)]
                    retvec[i] = true
                    break
                end
            end
            if retvec[i] == true
                break
            end
        end
    end
    return retvec
end

"""
    identify_condition2(df::AbstractDataFrame,cols::Vector{Symbol},codes::Vector{String})

produces a vector of Boolean `true` or `false` depending on whether columns `cols` contain
any of the string values in `codes`. This function is designed to easily identify comorbid
conditions in inpatient or outpatient records that are coded in ICD-9 or ICD-10 codes.
The values in `cols` are compared to the values in `codes` exactly as they are specified.
"""
function identify_condition2(ip::AbstractDataFrame,vicd::Vector{Symbol},codes::Vector{String})

    retvec = zeros(Bool,nrow(ip))

    for i in 1:nrow(ip)
        for v in vicd
            if ismissing(ip[i,v])
                continue
            end

            if ip[i,v] in codes
                retvec[i] = true
            end
        end
    end
    return retvec
end
