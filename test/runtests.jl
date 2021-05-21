module SpeckleTest

using Speckles
using Test

################################################################################
function TestSpeckles()
    @testset "Speckles.jl" begin
        # Write your tests here.
    end
end

export TestSpeckles

################################################################################
function TestLightSource()
    @testset "LightSource.jl" begin
        include("testLightSource.jl")
    end
end

export TestLightSource

################################################################################
function TestDetector()
    @testset "Detector.jl" begin
        include("testDetector.jl")        
    end
end

export TestDetector

################################################################################
function TestAll()
    TestSpeckles()
    TestLightSource()
    TestDetector()
end

export TestAll

end

import .SpeckleTest

SpeckleTest.TestDetector()
