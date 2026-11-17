using LinearAlgebra
"""
Approximation de la solution du problème 

    min qₖ(s) = s'gₖ + 1/2 s' Hₖ s, sous la contrainte ‖s‖ ≤ Δₖ

# Syntaxe

    s, flag_gct = gct(g, H, Δ; kwargs...)

# Entrées

    - g : (Vector{<:Real}) le vecteur gₖ
    - H : (Matrix{<:Real}) la matrice Hₖ
    - Δ : (Real) le scalaire Δₖ
    - kwargs  : les options sous formes d'arguments "keywords", c'est-à-dire des arguments nommés
        • max_iter : le nombre maximal d'iterations (optionnel, par défaut 100)
        • tol_abs  : la tolérance absolue (optionnel, par défaut 1e-10)
        • tol_rel  : la tolérance relative (optionnel, par défaut 1e-8)

# Sorties

    - s : (Vector{<:Real}) une approximation de la solution du problème
    - flag_gct : (Int) indique la façon dont l'algorithme s'est arrêté
        • 0 : convergence à l'intérieur de la région de confiance (‖g_j‖ assez petit, y compris g = 0)
        • 1 : courbure négative (κ_j ≤ 0), le pas est sur le bord de la région de confiance
        • 2 : sortie de la région de confiance, le pas est sur le bord de la région de confiance
        • 3 : nombre maximal d'itérations atteint

# Exemple d'appel

    g = [0; 0]
    H = [7 0 ; 0 2]
    Δ = 1
    s, flag_gct = gct(g, H, Δ, max_iter = 10)

"""
function gct(g::Vector{<:Real}, H::Matrix{<:Real}, Δ::Real; 
    max_iter::Int = 100, 
    tol_abs::Real = 1e-10, 
    tol_rel::Real = 1e-8)

    s = zeros(length(g))
    flag_gct = -1

    # à compléter
   return s, flag_gct
   
end