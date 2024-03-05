using LevenbergMarquardtFit
using Test
using Aqua
using JET

@testset "LevenbergMarquardtFit.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(LevenbergMarquardtFit; ambiguities = false,)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(LevenbergMarquardtFit; target_defined_modules = true)
    end
    # Write your tests here.
end
