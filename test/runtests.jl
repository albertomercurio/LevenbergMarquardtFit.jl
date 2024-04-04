using LevenbergMarquardtFit
using Test
using Aqua
using JET
using LinearAlgebra

@testset "LevenbergMarquardtFit.jl" begin
    function f(x,p)
        @. p[1] * exp(-(x - p[3])^2 / $(2 * p[2]^2))
    end
    
    xdata = range(0, 10, length=10000) |> collect
    
    p0 = [1, 9, 0.1]
    p_correct = [2, 1, 5]
    ydata = f(xdata, (p_correct)) .+ 0.0001*randn(length(xdata))
    
    lower = [0, 0, 0.0]
    upper = [10.0, 20, 10.0]

    p_opt = curve_fit(f, xdata, ydata, p0, lower=lower, upper=upper).params

    ADM = AutoDiffMethod(f, p0, xdata);
    p_opt_auto = curve_fit(f, xdata, ydata, p0, lower=lower, upper=upper, diff_method=ADM).params

    @test p_opt ≈ p_correct atol=1e-5
    @test p_opt_auto ≈ p_correct atol=1e-5
end

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(LevenbergMarquardtFit; ambiguities = false,)
end

@testset "Code quality (JET.jl)" begin
    if VERSION >= v"1.9"
        JET.test_package(LevenbergMarquardtFit; target_defined_modules=true, ignore_missing_comparison=true)
    end
end