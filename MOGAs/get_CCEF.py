import numpy as np
from models.card_mean_var import solve_MV
from MOGAs.nsga2.nsga2 import *
from MOGAs.nsga2_rule_guided.nsga2_rule_guided import *
from result_analysis.plot_UEF import *

def compute_CCEF(assets_return, k, T, indtrack, nGenerations = 1000):

	#'''
	#the following dictionaries contain information about the results obtained in each run of the metaheuristic:

	time_solutions = dict([]) # dictionary containing CPUtime obtained in each run of the approx algorithm
	time_to_target_sol = dict([]) # dictionary containing time to obtain the target solution in each run of the approx algorithm
	f_solutions = dict([]) # dictionary containing f_val obtained in each run of the approx algorithm
	it_solutions = dict([]) # dictionary containing niters obtained in each run of the approx algorithm
	w_solutions = dict([]) # dictionary containing  w vector obtained in each run of the approx algorithm

	#- where: dict.key = "execution index"
	#'''

	r = np.array(assets_return)
	r_filtered_by_T = r[:,0:T] # consider a return matrix with T return samples

	# get expected return
	assets_expected_return = np.mean(r_filtered_by_T, axis = 1)
	assets_std_dev = np.std(r_filtered_by_T, axis = 1)

	# get corr coef between each asset pair
	assets_corr_mat = np.corrcoef(r_filtered_by_T)
	assets_var_mat = np.cov(r_filtered_by_T)

	# get unconstrained efficient frontier
	CCEF = nsga2(assets_expected_return, assets_var_mat, k, indtrack, nGenerations)

	#print("get cardinality constrained Eff. Frontier using NSGA-II")
	#time_solutions, f_solutions, it_solutions, time_to_target_sol,  w_solutions = GA_torrubiano(model, r, k, T, num_executions, target_sol)

	return CCEF

def compute_CCEF_rule_guided(filehandle_init, rule, gamma, assets_return, k, T, indtrack):

	#'''
	#the following dictionaries contain information about the results obtained in each run of the metaheuristic:

	time_solutions = dict([]) # dictionary containing CPUtime obtained in each run of the approx algorithm
	time_to_target_sol = dict([]) # dictionary containing time to obtain the target solution in each run of the approx algorithm
	f_solutions = dict([]) # dictionary containing f_val obtained in each run of the approx algorithm
	it_solutions = dict([]) # dictionary containing niters obtained in each run of the approx algorithm
	w_solutions = dict([]) # dictionary containing  w vector obtained in each run of the approx algorithm

	#- where: dict.key = "execution index"
	#'''

	r = np.array(assets_return)
	r_filtered_by_T = r[:,0:T] # consider a return matrix with T return samples

	# get expected return
	assets_expected_return = np.mean(r_filtered_by_T, axis = 1)
	assets_std_dev = np.std(r_filtered_by_T, axis = 1)

	# get corr coef between each asset pair
	assets_corr_mat = np.corrcoef(r_filtered_by_T)
	assets_var_mat = np.cov(r_filtered_by_T)

	print("get cardinality constrained Eff. Frontier using NSGA-II guided by DRSA rules")
	X = get_initial_sol(filehandle_init)
	gen = 1
	CCEF, X, feasibility = nsga2_rule_guided(gamma, X, rule, assets_expected_return, assets_var_mat, k, indtrack)
	return CCEF, X

def plot_frontier(CCEF, indtrack, k):
	frontier_data = dict([]) # col 1: x(var), col2: y(exp), col3: uef or nsga-ii
	var = []
	exp = []
	front_type = []
	for ind in range(len(CCEF)):
		var.append(CCEF[ind][1])
		exp.append(CCEF[ind][0])
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

def get_initial_sol(filehandle):
	# get sample of solutions for this execution
	ind_file = ''
	for line in filehandle:
		ind_file += line
	ind_file2 = ind_file.split(']') 
	R = []
	count = 0
	for ind_raw in ind_file2:
		ind_cleaned = ind_raw.split(',')
		R.append([float(ind) for ind in ind_cleaned])
		count += 1	
		
	return np.array(R)
