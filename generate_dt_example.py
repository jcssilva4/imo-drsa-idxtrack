
from data.DBInstances import Instance
from data_tables import *
#from get_solution import get_solution
from MOGAs.get_portfolio import get_solution_approx
#from result_analysis.read_cplex_result_file import get_f_val
from MOGAs.get_CCEF import *
from rulegen_test import *
import json
# plots
import matplotlib.pyplot as plt 
import seaborn as sns

# launch gateway server that will communicate with jRS lib
gateway = launch_gatewayserver()

# experiment 1

# model and ga parameters
instances = [6]
k = 0 # select cardinalities
card_inst = dict([]) # according to guastaroba2012
card_inst[1] = 10
card_inst[2] = 10
card_inst[3] = 10
card_inst[4] = 10
card_inst[5] = 10
card_inst[6] = 40
oosT = 190 ## out-of-sample returns sample size
T = 290 - oosT
# adjust oosT
oosT = 100 #190

# simulated IMO-DRSA paramaters
intrLvl_max = 4 # maximum interaction level (lvl1: init sol, lvl2: sol after first interaction,....)
DM_types = ['TE','ER','EQ'] # 3 options: 1- TE (prefers tracking error), 2- ER (preferes excess return), 3 - EQ (Equilibrated investor - equal weights) 
DM_weights = dict([])
DM_weights['TE'] = [0.6, 0.3, 0.05, 0.05]
DM_weights['ER'] = [0.05, 0.05, 0.6, 0.3]
DM_weights['EQ'] = [0.25, 0.25, 0.25, 0.25]
frontier_filter = 'closer'
DT_presentations = ['visual','par_quant','nonpar_quant'] # par_quant (parametric quantile), nonpar_quant(non-parametric quantile), visual
DT_size = 3 # Num. of decision examples per data table = 2*DT_size
# GA optimization time
max_time_S0 = 5 # max_time for nsga2
max_time_S = 20 # max_time for rule-guided nsga2
# GA rebalance configurations
max_time_rebal = max_time_S # max time for portfolio rebalance
budget_time_rebal = 0 #120 (after single period evaluation we will select the best frontier filter and DT_presentation)
max_num_of_rebals = int(budget_time_rebal/max_time_rebal)
# GA population: num of individuals 
nIndividuals = 100
exc = 1
for instNum in instances:
	# 1. Read the data 
	thisInstance = Instance(instNum, 'beasley')
	rets = np.array(thisInstance.returns_i) # assets returns 
	I = np.array(thisInstance.returns_Index) # index returns 
	# set the cardinality
	k = card_inst[instNum]
	# Each folder instance/DM_type/frontier_filter/exc/ contains two files: results_optim.csv and results_oos.csv
	# Each folder instance/DM_type/ contains one file: results_optim_nsga2.csv (in the nsga2 experiments, after comparing the best filtering technique, we can use the CPUtime of each IMO-DRSA interaction as an input to limit the optim time budget of nsga2)
	
	# start the simulated interaction and save the results in results_optim.csv
	rets_is = rets[:,0:T] # assets returns in-sample
	I_is = I[0:T] # assets returns in-sample
	results_optim = dict([])
	# start the simulation
	R_current = []
	selected_rule = []
	intrLvl = 1
	print("generating S" + str(intrLvl-1))
	rules_not_induced = True
	CCEF, R_current, CPUTime = compute_CCEF(nIndividuals, rets_is, I_is, k, max_time_S0)
	total_CPUTime = CPUTime
	# sample 2*DT_size solutions from the non-dominated front
	rank_idx, fit, trash = sample_rand_sols(rets_is, I_is, DT_size, R_current,'visual')
	# get the evaluation of each portfolio in each objective
	final_fit = dict([])
	final_fit['obj1'] = []
	final_fit['obj2'] = []
	final_fit['obj3'] = []
	final_fit['obj4'] = []
	final_fit['f(TE)'] = []
	final_fit['f(ER)'] = []
	final_fit['f(EQ)'] = []
	for selected_ind in rank_idx:
		fit_ind = fit[selected_ind]
		for fit_m in range(len(fit_ind)):
			final_fit['obj'+str(fit_m+1)].append(fit_ind[fit_m])
	examples = [] # evaluations of portfolios that compose this data table
	R_all = [] # used to get the best portfolio (choice)
	R_current = np.array(R_current)
	for idx in rank_idx:
		examples.append(R_current[idx,:].tolist())
		R_all.append(R_current[idx])
	for DM_type in DM_types:
		# now we need to calculate the pref_func of the examples
		pref_order, DT_examples, pref_func = get_pref_func(R_all, DM_weights[DM_type], rets_is, I_is)
		final_fit['f(' + DM_type + ')'] = [pfunc for pfunc in pref_func]
		for DT_presentation in DT_presentations:
			# set the file path and create it if it doesn't exist
			filepath = "results/examples_interaction/" + DT_presentation + "/"
			if not os.path.exists(filepath):
				os.makedirs(filepath) # create this path if it not exists
			print("writing: " + filepath)
			print("DT_presentation (in-sample): " + str(DT_presentation))
			samples = get_drsa_eval(np.array(R_all), rets_is, I_is, DT_presentation)
			samples = np.array(samples)
			bestAlternative, DT = write_this_table_examples(R_all, samples, DT_size, DM_weights[DM_type], DT_presentation, rets_is, I_is,filepath + DM_type)
			# induce a set of rules and sample a unique rule
			rules = get_rules(gateway)
			f = open(filepath + DM_type + "rules.txt", 'w')  
			f.write(str(rules))
			f.close()
			# save outputs for this DT_presentation
			if DT_presentation in ['visual']:
				# save png that was genetared with a pandas df
				results_sim = dict([])
				results_sim["t"] = []
				results_sim["cumRet"] = []
				results_sim["class"] = []
				results_sim["portfolio"] = []
				sol_id = 0
				for portfolio in R_all:
					cumRet = 1
					for t in range(rets_is.shape[1]):
						results_sim["portfolio"].append("x" + str(sol_id+1))
						results_sim["t"].append(t+1)
						final_class = 'other'
						if sol_id in pref_order[0:DT_size]:
							final_class = 'good'
						results_sim["class"].append(final_class)
						cumRet = cumRet*(1+(portfolio @ rets_is[:,t]))
						results_sim["cumRet"].append(cumRet)
					sol_id += 1
				cumRet = 1
				for t in range(I_is.shape[0]):
					results_sim["portfolio"].append('index')
					results_sim["t"].append(t+1)
					results_sim["class"].append('index')
					cumRet = cumRet*(1+I_is[t])
					results_sim["cumRet"].append(cumRet)
				finalsim_df = pd.DataFrame.from_dict(results_sim)
				sns_plot = sns.lineplot(data = finalsim_df, x="t", y="cumRet", hue = "class", units = 'portfolio', estimator = None, hue_order = ['good','other','index'])
				sns_plot.get_figure().savefig(filepath + DM_type + "cumRet.png")
				plt.close('all')

	filepath = "results/examples_interaction/"
	fit_df = pd.DataFrame.from_dict(final_fit)
	f = open(filepath + "final_fit.txt", 'w')  
	f.write(fit_df.to_latex(index=False))

								














