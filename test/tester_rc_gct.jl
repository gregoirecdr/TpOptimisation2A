using Test
include("../test/fonctions_de_tests.jl")

"""
Tester l'algorithme des régions de confiance avec le gradient conjugué tronqué

# Entrées

	- afficher : (Bool) affichage ou non des résultats de chaque test

# les cas de test (dans l'ordre)

	- fct 1 : x011 (1 itération : f1 est quadratique et la région de confiance contient le pas de Newton), x012
	- fct 2 : x021, x022, x023
	- priorité des flags quand le critère d'arrêt est atteint à la dernière itération autorisée
	- pas nul (max_iter_gct = 0) : pas rejeté et rayon réduit

"""
function tester_rc_gct(rc::Function, afficher::Bool)

    println("Affichage des résultats des algorithmes : ", afficher, "\n")

	# La tolérance utilisée dans les tests
	tol_erreur = 1e-4

	# 
	algo_pas = "gct"

	# l'ensemble de tests
	Test.@testset "RC et gct" begin

		# cas de test 1
		x_sol, f_sol, flag, nb_iters, _ = rc(fct1, grad_fct1, hess_fct1, pts1.x011, algo_pas=algo_pas)
		if (afficher)
			afficher_resultats("RC et " * algo_pas, "f1", pts1.x011, x_sol, f_sol, flag, nb_iters, sol_exacte_fct1)
		end
		Test.@test x_sol ≈ sol_exacte_fct1 atol=tol_erreur
		Test.@test flag == 0
		Test.@test nb_iters == 1

		# cas de test 2
		x_sol, f_sol, flag, nb_iters  = rc(fct1, grad_fct1, hess_fct1, pts1.x012, algo_pas=algo_pas)
		if (afficher)
			afficher_resultats("RC et " * algo_pas, "f1", pts1.x012, x_sol, f_sol, flag, nb_iters, sol_exacte_fct1)
		end
		Test.@test x_sol ≈ sol_exacte_fct1 atol=tol_erreur
		Test.@test flag == 0
		Test.@test nb_iters ≤ 6 # 3 pour l'implémentation de référence

		# cas de test 3
		x_sol, f_sol, flag, nb_iters, _ = rc(fct2, grad_fct2, hess_fct2, pts1.x021, algo_pas=algo_pas)
		if (afficher)
			afficher_resultats("RC et " * algo_pas, "f2", pts1.x021, x_sol, f_sol, flag, nb_iters, sol_exacte_fct2)
		end
		Test.@test x_sol ≈ sol_exacte_fct2 atol=tol_erreur
		Test.@test flag in (0, 1, 2) # convergence ou stagnation
		Test.@test nb_iters ≤ 62 # 31 pour l'implémentation de référence

		# cas de test 4
		x_sol, f_sol, flag, nb_iters, _ = rc(fct2, grad_fct2, hess_fct2, pts1.x022, algo_pas=algo_pas)
		if (afficher)
			afficher_resultats("RC et " * algo_pas, "f2", pts1.x022, x_sol, f_sol, flag, nb_iters, sol_exacte_fct2)
		end
		Test.@test x_sol ≈ sol_exacte_fct2 atol=tol_erreur
		Test.@test flag in (0, 1, 2) # convergence ou stagnation
		Test.@test nb_iters ≤ 90 # 44 pour l'implémentation de référence

		# cas de test 5
		x_sol, f_sol, flag, nb_iters, _ = rc(fct2, grad_fct2, hess_fct2, pts1.x023, algo_pas=algo_pas)
		if (afficher)
			afficher_resultats("RC et " * algo_pas, "f2", pts1.x023, x_sol, f_sol, flag, nb_iters, sol_exacte_fct2)
		end
		Test.@test x_sol ≈ sol_exacte_fct2 atol=tol_erreur
		Test.@test flag in (0, 1, 2) # convergence ou stagnation
		Test.@test nb_iters ≤ 40 # 19 pour l'implémentation de référence

		# cas de test 7 : les sorties supplémentaires Δs (rayons), flags_gct (type de sortie du gct) et ρs (rapports)
		# f1 depuis x011 : une seule itération, le rayon initial contient le pas de Newton et le gct converge à l'intérieur
		x_sol, f_sol, flag, nb_iters, xs, Δs, flags_gct = rc(fct1, grad_fct1, hess_fct1, pts1.x011, algo_pas=algo_pas)
		Test.@test length(Δs) == nb_iters + 1 == length(xs)
		Test.@test Δs[1] == 2 # Δ0 par défaut
		Test.@test flags_gct == [0]
		# f2 depuis x021 : un flag par itération, les rayons restent dans ]0, Δmax]
		x_sol, f_sol, flag, nb_iters, xs, Δs, flags_gct, ρs = rc(fct2, grad_fct2, hess_fct2, pts1.x021, algo_pas=algo_pas)
		Test.@test length(Δs) == nb_iters + 1 == length(xs)
		Test.@test length(flags_gct) == nb_iters
		Test.@test all(0 .< Δs .≤ 10)
		Test.@test all(flags_gct .∈ Ref(0:3))
		Test.@test flags_gct[end] == 0 # la dernière itération converge à l'intérieur de la région (pas de Newton)
		# ρs : un rapport par itération ; le pas est accepté (l'itéré change) si et seulement si ρk ≥ η1 = 0.25
		Test.@test length(ρs) == nb_iters
		Test.@test (ρs .≥ 0.25) == [xs[k+1] != xs[k] for k in 1:nb_iters]
		# avec une seule itération du gct, une partie des pas est tronquée par max_iter_gct
		_, _, _, _, _, _, flags_gct = rc(fct2, grad_fct2, hess_fct2, pts1.x021, algo_pas=algo_pas, max_iter_gct=1, max_iter=100)
		Test.@test all(flags_gct .∈ Ref(0:3))
		Test.@test 3 in flags_gct

		# cas de test 8 : décroissance prédite nulle (max_iter_gct = 0 donne un pas nul) : le pas est rejeté et le rayon réduit
		# (et non ρ = NaN, pour lequel ni l'itéré ni le rayon n'étaient modifiés)
		_, _, flag, nb_iters, xs, Δs, _, ρs = rc(fct1, grad_fct1, hess_fct1, pts1.x011, algo_pas=algo_pas, max_iter_gct=0, max_iter=10)
		Test.@testset "pas nul" begin
			Test.@test flag == 3
			Test.@test !any(isnan, ρs)
			Test.@test all(x -> x == pts1.x011, xs)
			Test.@test Δs[end] < Δs[1]
		end

		# cas de test 6 : le critère d'arrêt est atteint à la dernière itération autorisée.
		# On relève d'abord le nombre d'itérations et le flag avec les paramètres par défaut, puis on relance
		# avec max_iter égal à ce nombre : le flag doit rester celui du critère le plus prioritaire (et non 3).
		for (f, g, H, x0, sol) in ((fct1, grad_fct1, hess_fct1, pts1.x011, sol_exacte_fct1),
		                           (fct1, grad_fct1, hess_fct1, pts1.x012, sol_exacte_fct1),
		                           (fct2, grad_fct2, hess_fct2, pts1.x021, sol_exacte_fct2))
			_, _, flag_ref, nb_ref, _ = rc(f, g, H, x0, algo_pas=algo_pas)
			x_sol, f_sol, flag, nb_iters, _ = rc(f, g, H, x0, algo_pas=algo_pas, max_iter=nb_ref)
			Test.@test x_sol ≈ sol atol=tol_erreur
			Test.@test nb_iters == nb_ref
			Test.@test flag == flag_ref
		end
	
	end
end
