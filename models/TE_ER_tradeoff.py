import numpy as np

def get_fitness_individuals(R, rets, I):
	# [f1 = mean(te_t), f2 = max(te_t), f3 = mean(er_t), f4 = max(-er_t)]
	fit = []
	for individual in R:
		fit.append(get_objs_individual(individual,rets,I))
	return fit

def get_objs_individual(individual,rets,I):
	weights = individual[:rets.shape[0]]
	portfolio_ret = weights @ rets
	ER = portfolio_ret - I
	TE = np.abs(ER)
	meanTE = np.mean(TE)
	maxTE = np.max(TE)
	meanER = np.mean(ER)
	maxNegER = np.max(-ER)
	return [meanTE, maxTE, meanER, maxNegER]

def get_pref_dir():
	pref_dir = [-1, -1, 1, -1]
	return pref_dir

def get_normalization_coefs():
	f_max = [1,1,1,1]
	f_min = [0,0,-1,-1]
	return f_max, f_min