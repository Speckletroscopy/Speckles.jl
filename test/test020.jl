module test020
using Speckles

function run1()
    Speckles.run_020()
end
export run1

end

import .test020

test020.run1()
