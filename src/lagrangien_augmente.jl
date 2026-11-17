using LinearAlgebra
include("../src/newton.jl")
include("../src/regions_de_confiance.jl")
"""

Approximation d'une solution au problème 

    min f(x), x ∈ Rⁿ, sous la contrainte c(x) = 0,

par l'algorithme du lagrangien augmenté.

# Syntaxe

    x_sol, f_sol, flag, nb_iters, xs, μs, λs = lagrangien_augmente(f, ∇f, ∇²f, c, ∇c, ∇²c, x0; kwargs...)

# Entrées

    - f      : (Function) la fonction à minimiser
    - ∇f     : (Function) le gradient de f
    - ∇²f    : (Function) la hessienne de f
    - c      : (Function) la contrainte, à valeur dans R (une seule contrainte d'égalité)
    - ∇c     : (Function) le gradient de la contrainte c
    - ∇²c    : (Function) la hessienne de la contrainte c
    - x0     : (Vector{<:Real}) itéré initial
    - kwargs : les options sous formes d'arguments "keywords"
        • max_iter  : (Int) le nombre maximal d'iterations (optionnel, par défaut 1000)
        • tol_abs   : (Real) la tolérance absolue (optionnel, par défaut 1e-6)
        • tol_rel   : (Real) la tolérance relative (optionnel, par défaut 1e-6)
        • λ0        : (Real) le multiplicateur de lagrange associé à c initial (optionnel, par défaut 2)
        • μ0        : (Real) le facteur initial de pénalité de la contrainte (optionnel, par défaut 10)
        • τ         : (Real) le facteur d'accroissement de μ (optionnel, par défaut 2)
        • algo_noc  : (String) l'algorithme sans contrainte à utiliser (optionnel, par défaut "rc-gct")
            * "newton"    : pour l'algorithme de Newton
            * "rc-cauchy" : pour les régions de confiance avec pas de Cauchy
            * "rc-gct"    : pour les régions de confiance avec gradient conjugué tronqué

# Sorties

    - x_sol    : (Vector{<:Real}) une approximation de la solution du problème
    - f_sol    : (Real) f(x_sol)
    - flag     : (Int) indique le critère sur lequel le programme s'est arrêté
        • 0 : convergence
        • 1 : nombre maximal d'itération dépassé
    - nb_iters : (Int) le nombre d'itérations faites par le programme
    - xs       : (Vector{Vector{<:Real}}) les itérés
    - μs       : (Vector{<:Real}) tableau des valeurs prises par μk au cours de l'exécution
    - λs       : (Vector{<:Real}) tableau des valeurs prises par λk au cours de l'exécution

# Exemple d'appel

    f(x)=100*(x[2]-x[1]^2)^2+(1-x[1])^2
    ∇f(x)=[-400*x[1]*(x[2]-x[1]^2)-2*(1-x[1]) ; 200*(x[2]-x[1]^2)]
    ∇²f(x)=[-400*(x[2]-3*x[1]^2)+2  -400*x[1];-400*x[1]  200]
    c(x) =  x[1]^2 + x[2]^2 - 1.5
    ∇c(x) = 2*x
    ∇²c(x) = [2 0; 0 2]
    x0 = [1; 0]
    x_sol, _ = lagrangien_augmente(f, ∇f, ∇²f, c, ∇c, ∇²c, x0, algo_noc="rc-gct")

"""
function lagrangien_augmente(f::Function, ∇f::Function, ∇²f::Function, 
        c::Function, ∇c::Function, ∇²c::Function, x0::Vector{<:Real}; 
        max_iter::Int=1000, tol_abs::Real=1e-6, tol_rel::Real=1e-6,
        λ0::Real=2, μ0::Real=10, τ::Real=2, algo_noc::String="rc-gct")

    #
    x_sol = x0
    f_sol = f(x_sol)
    flag  = -1
    nb_iters = 0
    xs = [Float64.(x0)] # vous pouvez faire push!(xs, xk) pour créer le vecteur des itérés (xs isa Vector{Vector{Float64}})
    μs = [Float64(μ0)] # vous pouvez faire push!(μs, μk) pour créer le vecteur des μk (μs isa Vector{Float64})
    λs = [Float64(λ0)] # vous pouvez faire push!(λs, λk) pour créer le vecteur des λk (λs isa Vector{Float64})

    # à compléter
    return x_sol, f_sol, flag, nb_iters, xs, μs, λs

end
