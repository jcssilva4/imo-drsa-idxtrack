from result_analysis.analysis_tools import *
import numpy as np
from data.DBInstances import Instance
from MOGAs.nsga2.tools import *
from MOGAs.nsga2.operators import *
from MOGAs.nsga2_rule_guided.tools_rule_guided import get_feasibility_individuals
import pandas as pd
# plot libs
import seaborn as sns
import matplotlib.pyplot as plt


intrLvl = 2 # select an interaction level
DM_type = 'risk-prone' # or risk-averse
numRuns = 30 # number of executions of the algorithm
instNum = 5
oosT = 50 ## out-of-sample period
T = 290 - oosT
DT_size = 3

# get mean number of feasible solutions for each portfolio in each type of choice file

# 1. Read the problem instance data
thisInstance = Instance(instNum, 'beasley')

# set infeasibility metrics
results_infeas = dict([])
results_infeas['type'] = []
results_infeas['t'] = []
results_infeas['infeasibility'] = []

# set feasibility metrics
results_feas = dict([])
results_feas['number of feasible solutions'] = []
results_feas['type'] = []
results_feas['t'] = []

'''
# comparison of the data table generation method

# read choice files for intrLvl, DM_type in each run
portfolios_closer, rules_closer = read_choices('closer', intrLvl, instNum, numRuns, DM_type) # an array of portfolio choices that are closer to the simulated ideal solutions
portfolios_farther, rules_farther = read_choices('farther', intrLvl, instNum, numRuns, DM_type) # an array of portfolio choices that are farther from the simulated ideal solutions
nAssets = portfolios_closer.shape[1]

for out_t in range(oosT): # loop over the out-of-sample period

	# get exp return and variance
	r = np.array(thisInstance.returns_i)
	r_filtered_by_T = r[:, out_t + 1 : T + out_t + 1] # consider a return matrix with T return samples
	#print(str(out_t + 1) + " - " + str(T + out_t + 1))
	#print(r_filtered_by_T.shape[1])
	# get expected return
	r_exp = np.mean(r_filtered_by_T, axis = 1)
	# get corr coef between each asset pair
	r_var = np.cov(r_filtered_by_T)

	# get fitness for each type of choice
	fit_closer = np.array(get_fitness_individuals(portfolios_closer, nAssets, r_exp, r_var))
	fit_farther = np.array(get_fitness_individuals(portfolios_farther, nAssets, r_exp, r_var))

	# get feasibility for each type of choice
	feasibility_closer = get_feasibility_individuals_all_excs(portfolios_closer, nAssets, rules_closer, fit_closer)
	feasibility_farther = get_feasibility_individuals_all_excs(portfolios_farther, nAssets, rules_closer, fit_farther) # rules are the same for closer and farther 
	len_data_table = 2*DT_size

	# get infeasible performance metrics
	for exc in range(numRuns):
		current_ind = exc*len_data_table
		# closer
		results_infeas['t'].append(out_t + 1)
		results_infeas['type'].append('closer')
		results_infeas['infeasibility'].append(np.sum([f for f in feasibility_closer[current_ind:current_ind+len_data_table] if f < 1]))
		# farther
		results_infeas['t'].append(out_t + 1)
		results_infeas['type'].append('farther')
		results_infeas['infeasibility'].append(np.sum([f for f in feasibility_farther[current_ind:current_ind+len_data_table] if f < 1]))

	# get feasible metrics
	for exc in range(numRuns):
		current_ind = exc*len_data_table
		# closer
		results_feas['t'].append(out_t + 1)
		results_feas['type'].append('closer')
		results_feas['number of feasible solutions'].append(len([f for f in feasibility_closer[current_ind:current_ind+len_data_table] if f == 1]))
		# farther
		results_feas['t'].append(out_t + 1)
		results_feas['type'].append('farther')
		results_feas['number of feasible solutions'].append(len([f for f in feasibility_farther[current_ind:current_ind+len_data_table] if f == 1]))

# prepare results
results_ready_infeas = pd.DataFrame.from_dict(results_infeas)
results_ready_feas = pd.DataFrame.from_dict(results_feas)
# plot results
plt.rcParams.update({'font.size': 15})
sns.lineplot(x="t", y="number of feasible solutions",hue="type",data=results_ready_feas, ci = 'sd')
plt.show()
sns.lineplot(x="t", y="infeasibility",hue="type",data=results_ready_infeas, ci = 'sd')
plt.show()

'''

#'''
# check impact of the level of interaction

# read choice files for intrLvl, DM_type in each run
portfolios_closer, rules_closer = read_choices('closer', 2, instNum, numRuns, DM_type) # an array of portfolio choices that are closer to the simulated ideal solutions
portfolios_farther, rules_farther = read_choices('closer', 3, instNum, numRuns, DM_type) # an array of portfolio choices that are farther from the simulated ideal solutions
nAssets = portfolios_closer.shape[1]

for out_t in range(oosT): # loop over the out-of-sample period

	# get exp return and variance
	r = np.array(thisInstance.returns_i)
	r_filtered_by_T = r[:, out_t + 1 : T + out_t + 1] # consider a return matrix with T return samples
	#print(str(out_t + 1) + " - " + str(T + out_t + 1))
	#print(r_filtered_by_T.shape[1])
	# get expected return
	r_exp = np.mean(r_filtered_by_T, axis = 1)
	# get corr coef between each asset pair
	r_var = np.cov(r_filtered_by_T)

	# get fitness for each type of choice
	fit_closer = np.array(get_fitness_individuals(portfolios_closer, nAssets, r_exp, r_var))
	fit_farther = np.array(get_fitness_individuals(portfolios_farther, nAssets, r_exp, r_var))

	# get feasibility for each type of choice
	feasibility_closer = get_feasibility_individuals_all_excs(portfolios_closer, nAssets, rules_closer, fit_closer)
	feasibility_farther = get_feasibility_individuals_all_excs(portfolios_farther, nAssets, rules_closer, fit_farther) # rules are the same for closer and farther 
	len_data_table = 2*DT_size
	# get infeasible performance metrics
	for exc in range(numRuns):
		current_ind = exc*len_data_table
		# closer
		results_infeas['t'].append(out_t + 1)
		results_infeas['type'].append('y = 2')
		results_infeas['infeasibility'].append(np.sum([f for f in feasibility_closer[current_ind:current_ind+len_data_table] if f < 1]))
		# farther
		results_infeas['t'].append(out_t + 1)
		results_infeas['type'].append('y = 3')
		results_infeas['infeasibility'].append(np.sum([f for f in feasibility_farther[current_ind:current_ind+len_data_table] if f < 1]))

	# get feasible metrics
	for exc in range(numRuns):
		current_ind = exc*len_data_table
		# closer
		results_feas['t'].append(out_t + 1)
		results_feas['type'].append('y = 2')
		results_feas['number of feasible solutions'].append(len([f for f in feasibility_closer[current_ind:current_ind+len_data_table] if f == 1]))
		# farther
		results_feas['t'].append(out_t + 1)
		results_feas['type'].append('y = 3')
		results_feas['number of feasible solutions'].append(len([f for f in feasibility_farther[current_ind:current_ind+len_data_table] if f == 1]))

# prepare results
results_ready_infeas = pd.DataFrame.from_dict(results_infeas)
results_ready_feas = pd.DataFrame.from_dict(results_feas)
# plot results
plt.rcParams.update({'font.size': 15})
sns.lineplot(x="t", y="number of feasible solutions",hue="type",data=results_ready_feas, ci = 'sd')
plt.show()
sns.lineplot(x="t", y="infeasibility",hue="type",data=results_ready_infeas, ci = 'sd')
plt.show()
#'''






