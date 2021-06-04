################################################################################
# String management for Speckles.jl
################################################################################


"""
    makeName(params::Dict; excpt::Dict = Dict())

Returns a string based on parameters in the simulation dictionary.
"""
function makeName(params::Dict; excpt::Dict = Dict())
    date = today()
    mm = length(month(date)) == 1 ? string(0,month(date)) : month(date)
    dd = length(day(date)) == 1 ? string(0,day(date)) : day(date)
    out = string(year(date),mm,dd,"_")
    for key in keys(params)
        keyname = key
        val = params[key]
        if isa(val,Array)
            val = length(val)
            keyname = string("len-",keyname)
        end
        out = string(out,keyname,"=",val,"_")
    end
    return out
end

export makeName

"""
    function paramTable(params::Dict)

Returns a markdown table with given parameters
"""
function paramTable(params::Dict,names::Dict = Dict())
    out = "| Name | Value |\n|:---:|:---:|\n"
    for (key,value) in params
        paramName = string(key)
        if key in keys(names)
            paramName = string(names[key])
        end
        out = string(out,"|$paramName|$value|\n")
    end
    return out

end

export paramTable