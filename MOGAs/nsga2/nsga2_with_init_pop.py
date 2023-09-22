import numpy as np
import itertools as it
from MOGAs.nsga2.tools import *
from MOGAs.nsga2.operators import *
import seaborn as sns
import matplotlib.pyplot as plt
#from result_analysis.plot_UEF import *
import pandas as pd
from models.TE_ER_tradeoff import *
import time

def nsga2_with_init_pop(nIndividuals, init_pop_X, rets, I, k, max_time):

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
	R = pop
	
	gen = 1
	t_start = time.time()
	#for gen in range(nGenerations):
	while time.time() - t_start < max_time:
		#print("generation: " + str(gen + 1))

		##  new population P_(t+1) ##

		# get fitness of R elements
		fit = get_fitness_individuals(R, rets, I)
		# non-dominated fronts over R (non-dominated sorting)
		F = non_dominated_sort(fit, pref_dir) # - f1_pref_dir = -1 because f1 is risk (minimize) and f2_pref_dir = 1 because f2 is return (maximize)
		# initialize the new population P_(t+1) 
		P, lastFrontIdx = set_new_pop(F, nIndividuals)  # try to get better performance by not calculating any front when the new population is complete
		#f_max, f_min = get_normalization_coefs(F, fit, lastFrontIdx, len(pref_dir))
		f_max, f_min = get_normalization_coefs()
		# get crowd_distance rank and front rank
		cdist, frank = front_dist_rank(f_max, f_min, F, fit, lastFrontIdx, len(pref_dir))
		# set the final new population
		P = set_new_pop(F, nIndividuals, cdist, lastFrontIdx, P) 

		## new offspring Q_(t+1) ##
		# initialize R_(t+1)
		R_new = np.zeros((nIndividuals, 2*nAssets))
		for ind in range(nIndividuals):
			R_new[ind,:] = R[P[ind],:]
		# binary tournament selection
		mating_pool = bin_tournament_selection(P, frank, cdist)
		# crossover
		Q = uniform_crossover_chang(mating_pool, R, nAssets, k, nIndividuals, p_c)
		#Q = uniform_crossover(mating_pool, R, nAssets, k, nIndividuals, p_c)
		# mutation
		Q = santanna_mutation(Q, nAssets, nIndividuals, 1, p_m)
		# R_(t+1) = union(P_(t+1),Q_(t+1))
		R = np.concatenate((R_new, Q))
		'''
		print(" -- -- constraint satisfaction -- --")
		for ind in range(R.shape[0]):
			print("K = " + str(sum(R[ind,nAssets:2*nAssets])) + " and sum(w_i) = " + str(sum(R[ind,0:nAssets])) )
		'''
		gen += 1
	
	#print("generation: " + str(gen))
	#print("time in secs: " + str(time.time()-t_start))

	# get the frontier of the last population
	CCEF = []
	final_sol_set = []
	R_new_idx = 0
	for ind in P:
		CCEF.append(fit[ind])
		final_sol_set.append(R_new[R_new_idx,:nAssets])
		R_new_idx += 1

	#print("plot")
	#print(CCEF)
	#plot_frontier(CCEF, indtrack, k)

	return CCEF, final_sol_set

def plot_frontier(CCEF, indtrack, k):
	frontier_data = dict([]) # col 1: x(var), col2: y(exp), col3: uef or nsga-ii
	var = []
	exp = []
	front_type = []
	for ind in range(len(CCEF)):
		var.append(CCEF[ind][0])
		exp.append(CCEF[ind][1])
		front_type.append("NSGA-II Frontier for K = " + str(k))	

	uef_exp, uef_var = get_UEF_data(indtrack) 
	for ind in range(len(uef_exp)):
		var.append(uef_var[ind])
		exp.append(uef_exp[ind])
		front_type.append("Unconstrained Eff. Frontier indtrack " + str(indtrack))

	frontier_data['var'] = var
	frontier_data['exp'] = exp
	frontier_data['type'] = front_type

	# get data
	frontiers = pd.DataFrame.from_dict(frontier_data)
	print(frontiers)
	sns.scatterplot(data=frontiers, x="var", y="exp", hue="type")

	plt.show()





