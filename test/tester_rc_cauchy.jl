using Test
include("../test/fonctions_de_tests.jl")

"""
Tester l'algorithme des régions de confiance avec le pas de Cauchy

# Entrées

   - afficher : (Bool) affichage ou non des résultats de chaque test

# les cas de test (dans l'ordre)

	- fct 1 : x011, x012
	- fct 2 : x021, x022, x023
	- priorité des flags quand le critère d'arrêt est atteint à la dernière itération autorisée
	- nom d'algorithme invalide (erreur même si x0 est déjà solution)

"""
function tester_rc_cauchy(rc::Function, afficher::Bool)

    println("Affichage des résultats des algorithmes : ", afficher, "\n")

	# Les tolérances utilisées dans les tests (le pas de Cauchy converge lentement sur f2)
	tol_erreur_f1 = 1e-4
	tol_erreur_f2 = 5e-2

	# 
	algo_pas = "cauchy"

	# l'ensemble de tests
	Test.@testset "RC et pas de Cauchy " begin

		# cas de test 1
		x_sol, f_sol, flag, nb_iters, _ = rc(fct1, grad_fct1, hess_fct1, pts1.x011, algo_pas=algo_pas)
		if (afficher)
			afficher_resultats("RC et " * algo_pas, "f1", pts1.x011, x_sol, f_sol, flag, nb_iters, sol_exacte_fct1)
		end
		Test.@test x_sol ≈ sol_exacte_fct1 atol=tol_erreur_f1
		Test.@test flag in (0, 2) # convergence ou stagnation de f
		Test.@test nb_iters ≤ 70 # 33 pour l'implémentation de référence

		# cas de test 2
		x_sol, f_sol, flag, nb_iters  = rc(fct1, grad_fct1, hess_fct1, pts1.x012, algo_pas=algo_pas)
		if (afficher)
			afficher_resultats("RC et " * algo_pas, "f1", pts1.x012, x_sol, f_sol, flag, nb_iters, sol_exacte_fct1)
		end
		Test.@test x_sol ≈ sol_exacte_fct1 atol=tol_erreur_f1
		Test.@test flag in (0, 2) # convergence ou stagnation de f
		Test.@test nb_iters ≤ 70 # 34 pour l'implémentation de référence

		# cas de test 3
		x_sol, f_sol, flag, nb_iters, _ = rc(fct2, grad_fct2, hess_fct2, pts1.x021, algo_pas=algo_pas)
		if (afficher)
			afficher_resultats("RC et " * algo_pas, "f2", pts1.x021, x_sol, f_sol, flag, nb_iters, sol_exacte_fct2)
		end
		# convergence très lente du pas de Cauchy sur Rosenbrock depuis x021 : le nombre maximal d'itérations
		# (5000 par défaut) est atteint pour l'implémentation de référence (flag 3), une convergence plus rapide est acceptée
		Test.@test x_sol ≈ sol_exacte_fct2 atol=tol_erreur_f2
		Test.@test flag in (0, 2, 3)
		Test.@test nb_iters ≤ 5000

		# cas de test 4
		x_sol, f_sol, flag, nb_iters, _ = rc(fct2, grad_fct2, hess_fct2, pts1.x022, algo_pas=algo_pas)
		if (afficher)
			afficher_resultats("RC et " * algo_pas, "f2", pts1.x022, x_sol, f_sol, flag, nb_iters, sol_exacte_fct2)
		end
		Test.@test x_sol ≈ sol_exacte_fct2 atol=tol_erreur_f2
		Test.@test flag in (0, 2)
		Test.@test nb_iters ≤ 2000 # 864 pour l'implémentation de référence

		# cas de test 5
		x_sol, f_sol, flag, nb_iters, _ = rc(fct2, grad_fct2, hess_fct2, pts1.x023, algo_pas=algo_pas)
		if (afficher)
			afficher_resultats("RC et " * algo_pas, "f2", pts1.x023, x_sol, f_sol, flag, nb_iters, sol_exacte_fct2)
		end
		Test.@test x_sol ≈ sol_exacte_fct2 atol=tol_erreur_f2
		Test.@test flag in (0, 2, 3)
		Test.@test nb_iters ≤ 5000

		# cas de test 7 : les sorties supplémentaires Δs (rayons), flags_gct (-1 avec le pas de Cauchy) et ρs (rapports)
		x_sol, f_sol, flag, nb_iters, xs, Δs, flags_gct, ρs = rc(fct1, grad_fct1, hess_fct1, pts1.x011, algo_pas=algo_pas)
		Test.@test length(Δs) == nb_iters + 1 == length(xs)
		Test.@test Δs[1] == 2 # Δ0 par défaut
		Test.@test length(flags_gct) == nb_iters
		Test.@test all(flags_gct .== -1)
		Test.@test all(0 .< Δs .≤ 10)
		# ρs : un rapport par itération ; le pas est accepté (l'itéré change) si et seulement si ρk ≥ η1 = 0.25
		Test.@test length(ρs) == nb_iters
		Test.@test (ρs .≥ 0.25) == [xs[k+1] != xs[k] for k in 1:nb_iters]

		# cas de test 8 : un nom d'algorithme invalide est signalé même si x0 est déjà solution (la boucle n'est pas exécutée)
		Test.@testset "algo_pas invalide" begin
			Test.@test_throws ErrorException rc(fct1, grad_fct1, hess_fct1, sol_exacte_fct1, algo_pas="toto")
		end

		# cas de test 6 : le critère d'arrêt est atteint à la dernière itération autorisée.
		# On relève d'abord le nombre d'itérations et le flag avec les paramètres par défaut, puis on relance
		# avec max_iter égal à ce nombre : le flag doit rester celui du critère le plus prioritaire (et non 3).
		for (f, g, H, x0, sol, tol_erreur) in ((fct1, grad_fct1, hess_fct1, pts1.x011, sol_exacte_fct1, tol_erreur_f1),
		                                       (fct2, grad_fct2, hess_fct2, pts1.x022, sol_exacte_fct2, tol_erreur_f2))
			_, _, flag_ref, nb_ref, _ = rc(f, g, H, x0, algo_pas=algo_pas)
			x_sol, f_sol, flag, nb_iters, _ = rc(f, g, H, x0, algo_pas=algo_pas, max_iter=nb_ref)
			Test.@test x_sol ≈ sol atol=tol_erreur
			Test.@test nb_iters == nb_ref
			Test.@test flag == flag_ref
		end
	
	end
end
