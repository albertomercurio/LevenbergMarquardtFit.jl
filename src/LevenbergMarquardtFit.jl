module LevenbergMarquardtFit

using LinearAlgebra
using ForwardDiff
using FiniteDiff

using FiniteDiff: JacobianCache
using ForwardDiff: JacobianConfig

export FitResult, AutoDiffMethod, FiniteDiffMethod, curve_fit

abstract type DiffMethod end

struct AutoDiffMethod{JCT,gFun} <: DiffMethod
    jacobian_cache::JCT
    g::gFun
end

struct FiniteDiffMethod{JCT,gFun} <: DiffMethod
    jacobian_cache::JCT
    g::gFun
end

function FiniteDiffMethod(f, p0, xdata)
    g = (dy, p) -> copyto!(dy, f(xdata, p))
    FiniteDiffMethod(JacobianCache(copy(p0), f(xdata, p0)), g)
end

function AutoDiffMethod(f, p0, xdata)
    g = p -> f(xdata, p)
    AutoDiffMethod(JacobianConfig(g, p0), g)
end

struct FitResult{PT<:AbstractVector, RT<:Real, JT<:AbstractArray, CT<:Bool, IT<:Integer}
    params::PT
    resid::RT
    jacobian::JT
    converged::CT
    iterations::IT
end

# write a show method in the text/plain format for the FitResult struct
function Base.show(io::IO, r::FitResult)
    println(io, "FitResult:")
    println(io, "  params: ", r.params)
    println(io, "  resid: ", r.resid)
    println(io, "  jacobian: ", r.jacobian)
    println(io, "  converged: ", r.converged)
    println(io, "  iterations: ", r.iterations)
end



_jacobian!(J, f::Function, x0, diff_method::FiniteDiffMethod) = FiniteDiff.finite_difference_jacobian!(J, f, x0, diff_method.jacobian_cache)

_jacobian!(J, f::Function, x0, diff_method::AutoDiffMethod) = ForwardDiff.jacobian!(J, f, x0, diff_method.jacobian_cache)



apply_f!(f, dy, p, xdata, diff_method::FiniteDiffMethod) = copyto!(dy, f(xdata, p))

apply_f!(f, dy, p, xdata, diff_method::AutoDiffMethod) = f(xdata, p)

function curve_fit(f::Function, xdata::AbstractArray, ydata::AbstractArray, p0::AbstractArray; 
    lower::AbstractArray=fill(-Inf, length(p0)), upper::AbstractArray=fill(Inf, length(p0)),
    maxiter::Int = 1000,
    λ::LT = 10,
    abstol::ET = 1e-8,
    reltol::ET = 1e-5,
    diff_method::MyMethod = FiniteDiffMethod(f, p0, xdata)) where {MyMethod<:DiffMethod,LT<:Real,ET<:Real}

    # check that the boundaries are consistent
    idxs = similar(p0, Bool)
    @. idxs = lower > upper
    any(idxs) && throw(ArgumentError("lower > upper for some parameters"))
    @. idxs = p0 < lower
    any(idxs) && throw(ArgumentError("p0 < lower for some parameters"))
    @. idxs = p0 > upper
    any(idxs) && throw(ArgumentError("p0 > upper for some parameters"))

    g = diff_method.g

    # Inizializzazione dei parametri
    n = length(p0)
    m = length(ydata)

    dy = copy(ydata)

    p = copy(p0)
    dp = similar(p)

    J = similar(dy, m, n)
    JJ = similar(p, n, n)

    cache_mat1 = similar(dy, n, n)
    cache_vec1 = similar(dy, n)
    cache_vec2 = similar(p, n)
    

    # Calcolo del gradiente e della matrice hessiana
    _jacobian!(J, g, p, diff_method)

    iter = 0
    converged = false

    # Iterazioni dell'algoritmo
    while !converged && iter <= maxiter
        residual = ydata - apply_f!(f, dy, p, xdata, diff_method)

        mul!(cache_vec1, J', residual)
        copyto!(cache_vec2, cache_vec1)
        mul!(cache_mat1, J', J)
        copyto!(JJ, cache_mat1)

        JJ[diagind(JJ)] .+= λ

        ldiv!(dp, factorize(JJ), cache_vec2)

        # Aggiornamento dei parametri
        p_new = p + dp

        # Controllo dei boundaries
        @. idxs = p_new < lower
        p_new[idxs] .= lower[idxs]
        @. idxs = p_new > upper
        p_new[idxs] .= upper[idxs]

        # Controllo del criterio di arresto basato sulla variazione della funzione obiettivo
        if norm(dp) < abstol || norm(dp) < reltol * norm(p)
            converged = true
        end

        # Aggiornamento dei parametri dell'algoritmo
        residual_new = ydata - apply_f!(f, dy, p_new, xdata, diff_method)
        if dot(residual_new, residual_new) < dot(residual, residual)
            λ /= 10
            p = p_new
            _jacobian!(J, g, p, diff_method)
        else
            λ *= 10
        end

        iter += 1
    end

    residual = ydata - apply_f!(f, dy, p, xdata, diff_method)
    return FitResult(p, dot(residual, residual), J, converged, iter)
end

end
