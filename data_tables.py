import numpy as np
from MOGAs.nsga2.tools import *
from MOGAs.nsga2.operators import *
from MOGAs.nsga2_rule_guided.tools_rule_guided import get_feasibility_individuals, get_drsa_eval
#from result_analysis.plot_UEF import *
import pandas as pd
import random
from models.TE_ER_tradeoff import *

'''
Test this decision table input for jMAF before implementing

**ATTRIBUTES

+ expRet: (continuous), gain
+ var: (continuous), cost
+ d: [other, good], gain, decision

decision: d

**PREFERENCES

expRet: gain
var: cost

**EXAMPLES

portfolio1: var1, ret1, other
portfolio2: var2, ret2, good
portfolio3: var3, ret3, other

**END
'''

def get_datatable(rets, I, rule, DT_size, R_final, intrLvl, DM_weights, frontier_filter, DT_presentation):

	pref_dir = get_pref_dir() # -> pref_dir_m = -1 to minimize and pref_dir_m = 1 to maximize
	# get fitness for all individuals of current R_final
	sample = get_drsa_eval(np.array(R_final), rets, I, DT_presentation)
	sample = np.array(sample)
	F = non_dominated_sort(sample, pref_dir) 
	#print("final_front: " + str(F))
	#print("|F[0]|_togenerate = " + str(len(F[0])))
	samples_bestFront = []
	samples_bestFront_idxs = []
	samples_bestFront_R = []
	if(intrLvl >  1): # if we already have at least one rule...
		# get fitness and feasibility of each individual of the population
		feasibility = get_feasibility_individuals(R_final, rule, sample)
		rank_feasible = np.argsort(feasibility).tolist()
		rank_feasible.reverse()
		#print("F[0]: " + str(F[0]))
		#print("R: ")
		#print(np.array(R_final))
		for s in F[0]:
			if feasibility[s] == 1:
				samples_bestFront.append(sample[s].tolist())
				samples_bestFront_idxs.append(s)
				samples_bestFront_R.append(R_final[s])
		# if there are not enough feasible solutions, accept nondominated solutions with the lowest infeasibility
		feas_idx = 0
		while len(samples_bestFront) < 2*DT_size:
			this_idx = rank_feasible[feas_idx]
			if this_idx not in samples_bestFront_idxs:
				samples_bestFront.append(sample[this_idx])
				samples_bestFront_idxs.append(this_idx)
				samples_bestFront_R.append(R_final[this_idx])
			feas_idx += 1
		#print("samples_bestFront_idxs: " + str(samples_bestFront_idxs))
		#print("samples_bestFront_R: ")
		#print(np.array(samples_bestFront_R))
		# get ideal sol and build the data table
		sample_fitness = np.array(samples_bestFront)
		# find the ideal solution considering the current population and rule
		ideal_sol = get_ideal_sol(rule, sample_fitness)	
		# find the distance between each solution contained in a feasible sol. of F[0] and the ideal solution
		dist = get_ideal_sol_dist(ideal_sol, sample_fitness)
		# generate a sample of solutions closer to the ideal solution
		rank_idx =  np.argsort(dist).tolist()
		#print("rank_idx: " + str(rank_idx))
		# adjust rank_idx
		if frontier_filter == 'farther': # if we want solutions farther from the ideal solution
			rank_idx.reverse()
	else:
		# get fitness of each individual of the population
		for s in F[0]:
			samples_bestFront.append(sample[s].tolist())
			samples_bestFront_R.append(R_final[s])
		sample_fitness = np.array(samples_bestFront)

		# generate a random sample of size = 2*DT_size from F[0]
		rank_idx = random.sample([idx for idx in range(len(sample_fitness))], 2*DT_size)

	# get data table
	bestAlternative, DT = write_this_table(rank_idx[:2*DT_size], samples_bestFront_R, sample_fitness, DT_size, DM_weights, DT_presentation, rets, I) 	

	return bestAlternative, DT

def write_this_table(rank_idx, R, sample, DT_size, DM_weights, DT_presentation, rets, I):

	examples = [] # evaluations of portfolios that compose this data table
	R_all = [] # used to get the best portfolio (choice)
	for idx in rank_idx:
		examples.append(sample[idx,:].tolist())
		R_all.append(R[idx])
	# now we need to calculate the pref_func of the examples
	pref_order, DT_examples, pref_func = get_pref_func(R_all, DM_weights, rets, I)
	#print(pref_order)
	# adjust the preference order of the examples according to the pref_func
	examples = [examples[pref_idx] for pref_idx in pref_order]
	#print(np.array(examples))
	# get the best alternative (choice)
	best_alternative = np.array(R_all[pref_order[0]])
	#print("bestAlternative: " + str(best_alternative))

	# open the file to save the current data table
	data_table = open("rulegen/data/currentDT.isf","w")
	# write the header
	data_table.writelines("**ATTRIBUTES\n\n")
	criteria, sense, trash = get_head_DT_presentation(DT_presentation)
	for obj in range(len(criteria)):
		data_table.writelines("+ " + criteria[obj] + ": (continuous)," + sense[obj] + "\n")
	data_table.writelines("+ d: [other, good], gain, decision\n")
	data_table.writelines("\ndecision: d\n")
	data_table.writelines("\n**PREFERENCES\n\n")
	for obj in range(len(criteria)):
		data_table.writelines(criteria[obj] + ": " + sense[obj] + "\n")
	data_table.writelines("\n**EXAMPLES\n")
	# write the decision examples
	exmplcounter = 1
	for obj in examples:
		#print(obj)
		line = '\nportfolio_' + str(exmplcounter) + ": "
		#print(line)
		for m in obj:
			#print(m)
			line += str(m)
			line += ', '
		if exmplcounter <= DT_size: # classify as good
			line += 'good'
		else:
			line += 'other'
		exmplcounter += 1
		#print(line)
		data_table.writelines(line)
	data_table.writelines("\n\n**END")

	return best_alternative, DT_examples

def get_pref_func(R, DM_weights, rets, I):
	pref_dir = get_pref_dir()
	examples = np.array(get_fitness_individuals(R, rets, I))
	#print(examples)
	pref_func = np.zeros(examples.shape[0])
	for m in range(examples.shape[1]):
		f_m_max = 1.0 #np.max(examples[:,m])
		if m in [0,1]: # TE, max(TE)
			f_m_min = 0.0 #np.min(examples[:,m])
		else: # ER, max(-ER)
			f_m_min = -1.0
		for exmpl in range(examples.shape[0]):
			pref_func[exmpl] += pref_dir[m]*DM_weights[m]*((examples[exmpl,m]-f_m_min)/(f_m_max-f_m_min))
	#print(pref_func)
	pref_order = np.argsort(pref_func).tolist()
	pref_order.reverse()
	return pref_order, examples.tolist(), pref_func

def get_pref_func_rebal(portfolio, DM_weights):
	pref_dir = get_pref_dir()
	pref_func = 0
	for m in range(portfolio.shape[0]):
		f_m_max = 1.0 #np.max(examples[:,m])
		if m in [0,1]: # TE, max(TE)
			f_m_min = 0.0 #np.min(examples[:,m])
		else: # ER, max(-ER)
			f_m_min = -1.0
		pref_func += pref_dir[m]*DM_weights[m]*((portfolio[m]-f_m_min)/(f_m_max-f_m_min))
	return pref_func

def get_head_DT_presentation(DT_presentation):
	criteria = []
	sense = []
	crit_dict = dict([])
	if DT_presentation == 'visual':
		criteria.append("meanTE")
		criteria.append("maxTE")
		criteria.append("meanER")
		criteria.append("maxNegER")
		crit_dict["meanTE"] = 0
		crit_dict["maxTE"] = 1
		crit_dict["meanER"] = 2
		crit_dict["maxNegER"] = 3
		sense.append("cost")
		sense.append("cost")
		sense.append("gain")
		sense.append("cost")
	if DT_presentation in ['par_quant','nonpar_quant']:
		criteria.append("TE1")
		criteria.append("TE50")
		criteria.append("ER1")
		criteria.append("ER50")
		criteria.append("ER99")
		crit_dict["TE1"] = 0
		crit_dict["TE50"] = 1
		crit_dict["ER1"] = 2
		crit_dict["ER50"] = 3
		crit_dict["ER99"] = 4
		sense.append("cost")
		sense.append("cost")
		sense.append("gain")
		sense.append("gain")
		sense.append("gain")
	return criteria, sense, crit_dict


# this function samples one rule from the set of rules that characterizes 'good' solutions
def sample_rule(rules, DT_presentation):
	trsh1, trsh2, crit_dict = get_head_DT_presentation(DT_presentation)
	# assemble rules that characterize 'good' solutions from the data table
	rule_set = [] # final rule set: list of dirs
	for rule in rules:
		_rule = dict([]) # three keys : ['obj', 'cond_dir', 'value']
		_rule['obj'] = []
		_rule['cond_dir'] = []
		_rule['value'] = []
		[conds_all, dec] = rule.split("=>")
		if 'good' in dec: # if this is a rule that classifies examples as 'good'
			conds = conds_all.split("&")
			for cond in conds:
				# get the associated objective
				for criterion in crit_dict.keys():
					if criterion in cond:
						_rule['obj'].append(crit_dict[criterion])
				# get the cond_dir
				cond_dir = '>='
				if "<=" in cond:
					cond_dir = '<='
				_rule['cond_dir'].append(cond_dir)
				# get the value
				val_raw = cond.split(cond_dir)
				val = val_raw[1].replace(")","")
				_rule['value'].append(float(val))
			#print(_rule)
			rule_set.append(_rule)
	#sample one rule from rule_set
	sample_rule = random.sample([rule for rule in rule_set], 1)
	return sample_rule[0]

def get_ideal_sol(rule, sample):

	rule_objs = rule['obj']
	rule_condDir = rule['cond_dir']
	# set the best examples for each conditional element of the rule
	ideal_sol = dict([])
	for cond_idx in range(len(rule_condDir)): # loop over all rule conditions (constraints of the problem)
		obj = rule_objs[cond_idx] # get the criterion associated with this restriction
		# get solution pref order with respect to this obj
		#print(sample[:,obj])
		idx_ = []
		idx_ = np.argsort(sample[:,obj]).tolist() # order this list in ascending order
		if rule_condDir[cond_idx] == '>=':
			idx_.reverse()
		
		ideal_sol[obj] = sample[idx_[0], obj]
	#print(ideal_sol)
	return ideal_sol

def get_ideal_sol_dist(ideal_sol, sample):

	sample_adapted = sample.copy()
	ideal_sol_ready = []
	dist = [] # distance between a feasible solution and the ideal solution
	for m in range(sample.shape[1]): # loop over all the objectives
		if m not in ideal_sol.keys():
			sample_adapted[:,m] = 0
			ideal_sol_ready.append(0)
		else:
			ideal_sol_ready.append(ideal_sol[m])
	# calculate dist from ideal solution
	for sol in sample_adapted:
		dist.append(np.linalg.norm(sol-ideal_sol_ready))
	#print("ideal_sol: " + str(ideal_sol_ready))
	#print("order: " + str(np.argsort(dist)))
	#print(sample_adapted)
	return dist

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

def get_bestAlternative_rebal(rets, I, rule, DT_size, R_final, frontier_filter, DT_presentation):

	pref_dir = get_pref_dir() # -> pref_dir_m = -1 to minimize and pref_dir_m = 1 to maximize
	# get fitness for all individuals of current R_final
	sample = get_drsa_eval(np.array(R_final), rets, I, DT_presentation)
	sample = np.array(sample)
	F = non_dominated_sort(sample, pref_dir) 
	#print("final_front: " + str(F))
	#print("|F[0]|_togenerate = " + str(len(F[0])))
	samples_bestFront = []
	samples_bestFront_idxs = []
	samples_bestFront_R = []
	# get fitness and feasibility of each individual of the population
	feasibility = get_feasibility_individuals(R_final, rule, sample)
	rank_feasible = np.argsort(feasibility).tolist()
	rank_feasible.reverse()
	for s in F[0]:
		if feasibility[s] == 1:
			samples_bestFront.append(sample[s].tolist())
			samples_bestFront_idxs.append(s)
			samples_bestFront_R.append(R_final[s])
	# if there are not enough feasible solutions, accept nondominated solutions with the lowest infeasibility
	feas_idx = 0
	while len(samples_bestFront) < 2*DT_size:
		this_idx = rank_feasible[feas_idx]
		if this_idx not in samples_bestFront_idxs:
			samples_bestFront.append(sample[this_idx])
			samples_bestFront_idxs.append(this_idx)
			samples_bestFront_R.append(R_final[this_idx])
		feas_idx += 1
	# get ideal sol and build the data table
	sample_fitness = np.array(samples_bestFront)
	# find the ideal solution considering the current population and rule
	ideal_sol = get_ideal_sol(rule, sample_fitness)	
	# find the distance between each solution contained in a feasible sol. of F[0] and the ideal solution
	dist = get_ideal_sol_dist(ideal_sol, sample_fitness)
	# generate a sample of solutions closer to the ideal solution
	rank_idx =  np.argsort(dist).tolist()
	# adjust rank_idx
	if frontier_filter == 'farther': # if we want solutions farther from the ideal solution
		rank_idx.reverse()

	# get best rebal alternative based on the distance w.r.t the ideal sol
	idx_best = rank_idx[0]
	bestAlternative = samples_bestFront_R[idx_best]

	return bestAlternative















