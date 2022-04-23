import numpy as np
from samples.data_tables import read_rule


def read_choices(_file, intrLvl, instNum, numRuns, DM_type):
	portfolios = np.array([])
	rules = [] # list of rules (or list of dicts, since each rule is a dict)
	for exc in range(1,numRuns+1):
		filepath = "samples/indtrack" + str(instNum) + "/" + str(DM_type) + "/exec" + str(exc) + "/"
		filepath_choice = filepath + "data_table" + str(intrLvl)  + "_choice_" + _file + ".txt"
		print("reading choices from " + filepath_choice)
		filehandle = open(filepath_choice, 'r')
		portfolio_choices = get_solutions(filehandle)
		if portfolios.shape[0] < 1: # initialize portfolios array
			portfolios = portfolio_choices
		else:
			portfolios = np.concatenate((portfolios, portfolio_choices), axis=0)

		# 2. Read the rule for this instNum, DM_type, exc, intrLvl (it is better to make a vector of rules...)
		rule = read_rule(instNum, DM_type, exc, intrLvl)
		for ind in range(portfolio_choices.shape[0]): # add 2*DM_type rules
			rules.append(rule)

	return portfolios, rules

def get_solutions(filehandle):
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


def get_feasibility_individuals_all_excs(R, nAssets, rules, fit):

	feasibility = []

	for ind_idx in range(R.shape[0]): # loop over all individuals
		rule_violation = 1 # assume that individual is a feasible solution (all feasible solutions have the same constr_violation = 1)
		rule = rules[ind_idx]
		#print(rule)
		rule_objs = rule['obj']
		rule_conds = rule['cond_dir']
		rule_vals = rule['value'] 
		for cond_idx in range(len(rule_conds)): # loop over all rule conditions (constraints of the problem)
			obj = 0 # assume the current obj is var
			rhs = rule_vals[cond_idx] # get right hand side of the constraint

			# compute the current constraint value for this individual
			if rule_objs[cond_idx] == 'expRet':
				obj = 1
			current_constr_val = fit[ind_idx][obj] 

			# calculate the constr violation for this individual in this rule condition
			constr_violation = 0
			if rule_conds[cond_idx] == '<=':
				constr_violation = rhs - current_constr_val
			else:
				constr_violation = current_constr_val - rhs

			# check feasibility for this condition
			if(constr_violation < 0): # infeasible solution 
				if rule_violation == 1: # initialize rule violation
					rule_violation = constr_violation 
				else: # increment rule violation
					rule_violation += constr_violation 

		feasibility.append(rule_violation)
	return feasibility


