module Test2

using IterTools
using Dates

function paramVector(params::Dict)
    k = collect(keys(params))
    v = collect(Iterators.product(collect(values(params))...))
    return map(vals->Dict(collect(zip(k,vals))),ivec(v))
end

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

function run()
    paramDict = Dict(
                     :tres=>[0.01,0.10],
                     :tmax=>[20.0],
                     :bigN=>[10,100],
                     :ωM=>[[456811.0,456815.0]],
                     :temp=>[1000],
                     :nbar=>[100],
                     :ntot=>[5000]
                    )
    paramVec = paramVector(paramDict)
    makeName(paramVec[1])
end

end

import .Test2
Test2.run()
