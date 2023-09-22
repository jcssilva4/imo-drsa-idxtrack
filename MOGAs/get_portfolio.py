import numpy as np

def get_solution_approx(assets_return, k, T, num_executions):

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

	#print("get cardinality constrained Eff. Frontier using NSGA-II")
	#time_solutions, f_solutions, it_solutions, time_to_target_sol,  w_solutions = GA_torrubiano(model, r, k, T, num_executions, target_sol)

	return time_solutions, f_solutions, it_solutions, time_to_target_sol, w_solutions