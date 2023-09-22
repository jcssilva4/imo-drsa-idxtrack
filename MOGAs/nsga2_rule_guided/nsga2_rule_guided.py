import numpy as np
import itertools as it
from MOGAs.nsga2_rule_guided.tools_rule_guided import *
from MOGAs.nsga2_rule_guided.operators import *
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from models.TE_ER_tradeoff import *
import time


def nsga2_rule_guided(nIndividuals, init_pop_X, rule, rets, I, k, DT_presentation, max_time):

	CCEF = dict([])
	ub,lb=1,0 

	# define your optimisation problem
	nAssets = rets.shape[0]
	pref_dir = get_pref_dir()

	# set real coded NSGA-II parameters
	p_c = 0.9 # crossover probability
	p_m = 1/nAssets # mutation probability


	# get initial population
	pop =  get_initial_pop_X(init_pop_X, k, nAssets)
	#pop = get_initial_pop_prefconstr(init_pop_X,nIndividuals, nAssets, k, r_exp, r_var, rule)
	R = pop
	gen = 1
	X_f = [] # X's feasible solutions
	F0_final = [] # X's non-dominated feasible solutions
	t_start = time.time()
	#while len(F0_final)/nIndividuals < gamma_current:
	while time.time() - t_start < max_time:
		##  new population P_(t+1) ##
		# get fitness and feasibility of R elements
		fit = get_fitness_individuals(R, rets, I)
		drsa_eval =  get_drsa_eval(R, rets, I, DT_presentation)
		feasibility = get_feasibility_individuals(R, rule, drsa_eval)
		# non-dominated fronts over R (non-dominated sorting)
		F = non_dominated_sort(fit, pref_dir) # - f1_pref_dir = -1 because f1 is risk (minimize) and f2_pref_dir = 1 because f2 is return (maximize)
		# initialize the new population P_(t+1) 
		P, lastFrontIdx = set_new_pop(F, nIndividuals)  # try to get better performance by not calculating any front when the new population is complete
		#f_max, f_min = get_normalization_coefs(F, fit, lastFrontIdx, len(pref_dir))
		f_max, f_min = get_normalization_coefs()
		# get crowd_distance rank and front rank
		cdist, frank = front_dist_rank(f_max, f_min, F, fit, lastFrontIdx, len(pref_dir))
		# set the final new population
		#P = set_new_pop(F, nIndividuals, cdist, lastFrontIdx, P) 
		P = set_new_pop(F, nIndividuals, feasibility, lastFrontIdx, P) 

		## new offspring Q_(t+1) ##

		# initialize R_(t+1)
		R_new = np.zeros((nIndividuals, 2*nAssets))
		feasibility_final = []
		fit_final = []# remove
		for ind in range(nIndividuals):
			feasibility_final.append(feasibility[P[ind]])
			R_new[ind,:] = R[P[ind],:]
			fit_final.append(fit[P[ind]])

		F0_final = []
		#for ind in F[0]:
		for ind in P:
			if feasibility[ind] == 1:
				F0_final.append(ind)


		# binary tournament selection
		mating_pool = bin_tournament_selection(P, frank, feasibility, cdist)
		# crossover
		Q = uniform_crossover_chang(mating_pool, R, nAssets, k, nIndividuals, p_c)
		#Q = uniform_crossover(mating_pool, R, nAssets, k, nIndividuals, p_c)
		# mutation
		Q = santanna_mutation(Q, nAssets, nIndividuals, 1, p_m)

		# R_(t+1) = union(P_(t+1),Q_(t+1))
		R = np.concatenate((R_new, Q))

		X_f = get_feasible_sols(feasibility_final)
		gen += 1

	#print("generation: " + str(gen))
	#print("time in secs: " + str(time.time()-t_start))
	#print("|X_f_F[0]|/|X|: " + str(len(F0_final)/nIndividuals))
	#print("feasible -> |F[0]| = " + str(len(F0_final)))

	# get the frontier of the last population
	CCEF = []
	final_feasible = []
	final_sol = []
	count_idx = 0
	for ind in P:
		CCEF.append(fit[ind])
		final_feasible.append(feasibility[ind])
		final_sol.append(R_new[count_idx,:nAssets])
		count_idx += 1
	#print("plot")
	#print(CCEF)
	#plot_frontier(CCEF, indtrack, k)
	return CCEF, np.array(final_sol), final_feasible








