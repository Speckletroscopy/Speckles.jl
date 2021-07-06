module MakeNB
using Weave

function run()
    filename = joinpath(pwd(),"simnb.jmd")
    nb_path  = joinpath(pwd(),"report")
    # Julia markdown to PDF
    weave(filename; doctype = "md2pdf", out_path = nb_path)
end
export run

end

import .MakeNB
MakeNB.run()

