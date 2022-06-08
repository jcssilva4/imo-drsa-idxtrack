import random
import numpy as np
from models.TE_ER_tradeoff import *
from scipy.stats import norm


def get_initial_pop_X(X, k, nAssets):

	initialPop = np.zeros((X.shape[0], 2*nAssets))

	for individual in range(X.shape[0]):

		temp_var = X[individual,:]

		final_sol = repair_vars(temp_var, k, nAssets)
		initialPop[individual, :] = final_sol[:]

	return initialPop
'''
def get_initial_pop_prefconstr(X, nIndividuals, nAssets, k, r_exp, r_var, rule):

	final_initial_pop = []

	initialPop = np.zeros((X.shape[0], 2*nAssets))

	for individual in range(X.shape[0]):

		temp_var = X[individual,:]

		final_sol = repair_vars(temp_var, k, nAssets)
		initialPop[individual, :] = final_sol[:]

	# check feasible individuals in this current population
	fit = get_fitness_individuals(initialPop, nAssets, r_exp, r_var)
	feasibility = get_feasibility_individuals(initialPop, nAssets, rule, fit)
	sorted_feasible = sorted( [[x,i] for (i,x) in enumerate(feasibility)], reverse=True )[:30] 
	print(sorted_feasible)

	while(len(final_initial_pop) < nIndividuals):

		# gather feasible individuals
		while len(final_initial_pop) < nIndividuals:
			for ind in sorted_feasible:
				final_initial_pop.append(initialPop.shape[0])
			individual += 1
		print("current init pop progress: " + str(len(final_initial_pop)/nIndividuals))

	return np.array(final_initial_pop)
'''
def repair_vars(temp_var, k, nAssets, crossover = False, bug_ver = False):

	if crossover:
		#  when crossover occurs, some solutions may have length != k (santanna2017)
		current_psize = np.sum(temp_var[nAssets:2*nAssets])
		diff = current_psize - k # current violation magnitude
		#print("diff: " + str(diff))
		if diff > 0: # remove assets that contain small weights
			temp_var = repair_vars(temp_var[:nAssets], k, nAssets, bug_ver = True)

		elif diff < 0: # add assets
			idx_add = random.sample([i for i in range(nAssets, 2*nAssets)],nAssets) 
			count = 0
			while diff != 0:
				if(temp_var[idx_add[count]] == 0): # include this
					temp_var[idx_add[count]] = 1
					diff += 1
				count += 1

		return temp_var

	repaired = []
	# get the K top values 
	sorted_vec = sorted( [[x,i] for (i,x) in enumerate(temp_var)], reverse=True )[:k] 
	kth_largest_weight = sorted_vec[k-1]
	#print(kth_largest_weight)

	# only the k highest weights will be maintained
	w_var = [] # continuous vars
	z_var = [] # binary vars
	k_constr = 0
	for var in range(nAssets):
		value = 0
		value_bin = 0

		if temp_var[var] >= kth_largest_weight[0] and k_constr < k:
			value = temp_var[var]
			value_bin = 1
			k_constr += 1
		w_var.append(value)
		z_var.append(value_bin)

	repaired = w_var
	repaired.extend(z_var)

	return repaired

def repair_weights(temp_var, nAssets):

	repaired = []
	weights = temp_var[0:nAssets]
	sum_weights = np.sum(weights)
	weights_normalized = np.dot((1/sum_weights),weights)
	repaired = weights_normalized.tolist()
	repaired.extend(temp_var[nAssets:])

	return repaired


def non_dominated_sort(fit, pref_dir):

	F = [] # non-dominated fronts
	evaluating = [e for e in range(len(fit))]

	while(len(evaluating)>0): # while there are elements being evaluated

		f_temp = [] # initialize a front
		for i in evaluating:
			# check if i is dominated by any element j in evaluating
			nonDominated = True # assume i is a non-dominated element
			j = 0
			while nonDominated and j < len(evaluating):
				eval_i = fit[i]
				eval_j = fit[evaluating[j]]
				count_dom = 0
				for obj in range(len(pref_dir)): # loop over all objectives
					if(pref_dir[obj]*eval_i[obj] < pref_dir[obj]*eval_j[obj]): # if i is dominated in this obj
						count_dom += 1
				if(count_dom == len(pref_dir)): # check if i is dominated by at least one j
					nonDominated = False
				j += 1
			if(nonDominated): # if i is not dominated
				f_temp.append(i) # add i to the current non-dominated front

		# remove elements that were included in the current front
		eval_temp = [] #remaining elements (dominated by elements of f_temp)
		for i in evaluating:
			if not i in f_temp:
				eval_temp.append(i)
		evaluating = eval_temp

		# update the non-dominated fronts vec
		F.append(f_temp)

		#
		# save some computation by stopping the algorithm earlier
		#new_P = [] # simulates a new pop
		#for f in F: # loop over all current fronts
		#	new_P.extend(f)
		#if(len(new_P) >= nIndividuals): 
		#	evaluating = 0 # the number of fronts are sufficient to compute the new population
		#

	return F


def get_feasibility_individuals(R, rule, drsa_eval):
	feasibility = []
	rule_objs = rule['obj']
	rule_conds = rule['cond_dir']
	rule_vals = rule['value']
	for ind_idx in range(R.shape[0]): # loop over all individuals
		rule_violation = 1 # assume that individual is a feasible solution (all feasible solutions have the same constr_violation = 1)
		for cond_idx in range(len(rule_conds)): # loop over all rule conditions (constraints of the problem)
			rhs = rule_vals[cond_idx] # get right hand side of the constraint
			obj = rule_objs[cond_idx] # get the associated objective
			current_constr_val = drsa_eval[ind_idx][obj] # get the current constr. value
			# calculate the constr violation for this individual in this rule condition
			constr_violation = 0
			if rule_conds[cond_idx] == '<=':
				constr_violation = rhs - current_constr_val
			else:
				constr_violation = current_constr_val - rhs
			# check feasibility for this condition
			if(constr_violation < 0): # infeasible solution 
				#print("infeasible: " + str(constr_violation))
				if rule_violation == 1: # initialize rule violation
					rule_violation = constr_violation 
				else: # increment rule violation
					rule_violation += constr_violation 
		feasibility.append(rule_violation)
	return feasibility

def get_feasible_sols(feasibility):
	feasible_sol_idxs = [] # stores indexes of feasible solutions 
	count = 0
	#print(feasibility)
	for f in feasibility:
		if f >= 0:
			#print("feasible")
			feasible_sol_idxs.append(count)
		count += 1
	return feasible_sol_idxs


def set_new_pop(F, nIndividuals, cdist = None, lastFrontIdx = None, P = []):

	if lastFrontIdx == None:
		#print("F: " + str(F))
		new_P = []
		i = 0
		while(len(new_P) + len(F[i]) < nIndividuals):
			new_P.extend(F[i])
			i += 1

		return new_P, i

	else:
		#print("P: " + str(P))
		#print("cdist: " + str(cdist))
		new_P = P
		remaining_elements = nIndividuals - len(P)
		cdist_unranked = []
		front = F[lastFrontIdx]
		for p in front:
			cdist_unranked.append(cdist[p])

		#print("cdist_unranked: " + str(cdist_unranked))

		# set the final new_P
		sorted_idx = [e[0] for e in sorted(enumerate(cdist_unranked), key=lambda x:x[1])]
		#print("sorted_idx: " + str(sorted_idx))
		for i in range(remaining_elements):
			new_P.append(front[sorted_idx[len(sorted_idx) - (i+1)]])

		#print("final_P: " + str(new_P))

		return new_P


'''
def get_normalization_coefs(F, fit, lastFrontIdx, nObjectives):

	temp_pop_idxs = []
	fit_obj = dict([])
	f_max = [] # contains the population-maximum for each obj
	f_min = [] # contains the population-minimum for each obj

	for i in range(lastFrontIdx + 1):
		temp_pop_idxs.extend(F[i])

	# initialize fit_obj
	for obj in range(nObjectives):
		fit_obj[obj] = []
	
	# get fit vectors for each obj func of solutions contained in F[0:lastFronIdx]
	for p in temp_pop_idxs:
		for obj in range(nObjectives):
			fit_obj[obj].append(fit[p][obj])

	# get f_max, f_min
	for obj in range(nObjectives):
		f_max.append(max(fit_obj[obj]))
		f_min.append(min(fit_obj[obj]))

	return [1,1], [0,0]

'''
def front_dist_rank(f_max, f_min, F, fit, lastFrontIdx, nObjectives):

	inf_const = 1e16 # represents a big number
	cdist = dict([]) # crowding distance
	frank = dict([]) # front rank
	for i in range(lastFrontIdx + 1): # loop over all fronts
		if(len(F[i]) < 3): # if F_i only have two solutions
			for p in F[i]:
				cdist[p] = inf_const 
				frank[p] = i + 1
		else: # F_i has intermediate solutions
			front = F[i]
			fit_obj = dict([])
			# initialize fit_obj
			for obj in range(nObjectives):
				fit_obj[obj] = []
			
			# index ranking for each obj 
			for p in front:
				frank[p] = i + 1
				cdist[p] = 0 # initialize cdist for the elements of this front
				for obj in range(nObjectives):
					fit_obj[obj].append(fit[p][obj])

			for obj in range(nObjectives):
				#sort front elements according to obj function m
				sorted_idx = [e[0] for e in sorted(enumerate(fit_obj[obj]), key=lambda x:x[1])]
				count = 0
				# get cdist for each element of this front
				for p in sorted_idx:
					if count == 0 or count == (len(front) - 1): # check if this is a border solution
						cdist[front[p]] = inf_const
					else:
						dist = (fit[front[sorted_idx[count+1]]][obj] - fit[front[sorted_idx[count-1]]][obj])/(f_max[obj] - f_min[obj])
						# get overall crowd dist
						if(dist > cdist[front[p]]):	
							cdist[front[p]] = dist
					count += 1

	return cdist, frank

def get_drsa_eval(R, rets, I, DT_presentation):
	drsa_eval = []
	prob1 = 0.99
	prob50 = 0.5
	prob99 = 0.01
	if DT_presentation == 'visual':
		drsa_eval = get_fitness_individuals(R, rets, I)
	if DT_presentation in ['par_quant','nonpar_quant']:
		for individual in R:
			weights = individual[:rets.shape[0]]
			portfolio_ret = weights @ rets
			ER = portfolio_ret - I
			TE = np.abs(ER)
			TE1 = 0
			TE50 = 0
			ER1 = 0
			ER50 = 0
			ER99 = 0
			minER = 0
			if DT_presentation == 'par_quant':
				TE1 = norm.ppf(prob1, loc=np.mean(TE), scale=np.std(TE))
				TE50 = norm.ppf(prob50, loc=np.mean(TE), scale=np.std(TE))
				ER1 = norm.ppf(prob1, loc=np.mean(ER), scale=np.std(ER))
				ER50 = norm.ppf(prob50, loc=np.mean(ER), scale=np.std(ER))
				ER99 = norm.ppf(prob99, loc=np.mean(ER), scale=np.std(ER))
			else:
				TE1 = np.quantile(TE,prob1)
				TE50 = np.quantile(TE,prob50)
				ER1 = np.quantile(ER,prob1)
				ER50 =  np.quantile(ER,prob50)
				ER99 = np.quantile(ER,prob99)
			drsa_eval.append([TE1, TE50, ER1, ER50, ER99])

	return drsa_eval
	

	