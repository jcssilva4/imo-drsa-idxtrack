

'''
Next:

1. check and test data_tables.py
'''

from data.DBInstances import Instance
from data_tables import *
#from get_solution import get_solution
from MOGAs.get_portfolio import get_solution_approx
#from result_analysis.read_cplex_result_file import get_f_val
from MOGAs.get_CCEF import *
from rulegen_test import *
import json

# launch gateway server that will communicate with jRS lib
gateway = launch_gatewayserver()

# experiment 1

# model and ga parameters
instances = [1,2,3,4,5,6]
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
DM_types = ['TE','ER','EQ'] # 3 options: 1- TE (prefers tracking error), 2- ER (prefers excess return), 3 - EQ (Equilibrated investor - equal weights)  
DM_weights = dict([])
DM_weights['TE'] = [0.6, 0.3, 0.05, 0.05]
DM_weights['ER'] = [0.05, 0.05, 0.6, 0.3]
DM_weights['EQ'] = [0.25, 0.25, 0.25, 0.25]
frontier_filters = ['closer','farther']
DT_presentations = ['visual','par_quant','nonpar_quant'] # par_quant (parametric quantile), nonpar_quant(non-parametric quantile), visual
numRuns = 30 # number of simulations per each DM type, instance, and frontier filter (each simulation consists in an interaction process from 1 to intrLvl_max).
DT_size = 3 # Num. of decision examples per data table = 2*DT_size
# GA optimization time
max_time_S0 = 10 # max_time for nsga2
max_time_S = 20 # max_time for rule-guided nsga2
# GA rebalance configurations
max_time_rebal = max_time_S # max time for portfolio rebalance
budget_time_rebal = 0 #120 (after single period evaluation we will select the best frontier filter and DT_presentation)
max_num_of_rebals = int(budget_time_rebal/max_time_rebal)
# GA population: num of individuals 
nIndividuals = 100
'''
#Experiment 1
for instNum in instances:
	# 1. Read the data 
	thisInstance = Instance(instNum, 'beasley')
	rets = np.array(thisInstance.returns_i) # assets returns 
	I = np.array(thisInstance.returns_Index) # index returns 
	# set the cardinality
	k = card_inst[instNum]
	for DM_type in DM_types:
		for frontier_filter in frontier_filters:
			for exc in range(1,numRuns+1): 

				# Each folder instance/DM_type/frontier_filter/exc/ contains two files: results_optim.csv and results_oos.csv
				# Each folder instance/DM_type/ contains one file: results_optim_nsga2.csv (in the nsga2 experiments, after comparing the best filtering technique, we can use the CPUtime of each IMO-DRSA interaction as an input to limit the optim time budget of nsga2)
				# set the file path and create it if it doesn't exist
				filepath = "results/core_experiments/indtrack" + str(instNum) + "/DM_" + DM_type + "/" + frontier_filter +"/sim" + str(exc) + "/"
				if not os.path.exists(filepath):
					os.makedirs(filepath) # create this path if it not exists
				print("writing: " + filepath)
				
				# start the simulated interaction and save the results in results_optim.csv
				rets_is = rets[:,0:T] # assets returns in-sample
				I_is = I[0:T] # assets returns in-sample
				results_optim = dict([])
				# simulation inputs
				results_optim["instNum"] = [] # 1,...,6
				results_optim["dm_type"] = [] # TE, ER or EQ
				results_optim["intrLvl"] = [] # interaction level: 1,2,...,intrLvl_max
				results_optim["frontier_filter"] = [] # closer or farther
				results_optim["simulation"] = [] # simulation number: exec
				results_optim["dt_type"] = [] # visual, par_ret, nonpar_ret
				# simulation outputs
				results_optim["current_population"] = [] # current population after the optimization (string)
				results_optim["optimTime"] = [] # time until satisfying the proportion of feasible individuals (gamma) in seconds
				results_optim["bestAlternative"] = [] # best alternative if the DM choose an alternative instead of choosing to continue the interaction
				results_optim["prefFunc"] = [] # evaluation of all the elemnts that constituted DT^y (list)
				results_optim["selected_rule"] = [] # rule selected in the current interaction (string)
				
				
				for DT_presentation in DT_presentations:
					print("DT_presentation (in-sample): " + str(DT_presentation))
					# start the simulation
					R_current = []
					selected_rule = []
					for intrLvl in range(1,intrLvl_max+1): # loop over all interaction levels
						print("generating S" + str(intrLvl-1))
						rules_not_induced = True
						total_CPUTime = 0
						while(rules_not_induced): # avoid cases in which no rules were induced or only rules that characterize 'other' types of portfolio
							if intrLvl == 1: # get the initial set of solutions S^0 + data_table1
								CCEF, R_current, CPUTime = compute_CCEF(nIndividuals, rets_is, I_is, k, max_time_S0)
							else: # generate S(intrLvl-1).txt using selected_rule and R_current
								CCEF, R_current, CPUTime = compute_CCEF_rule_guided(nIndividuals, R_current, selected_rule, rets_is, I_is, k,  DT_presentation, max_time_S)
							total_CPUTime += CPUTime
							# generate data_table(intrLvl)
							bestAlternative, DT, unique_sol_found = get_datatable(rets_is, I_is, selected_rule, DT_size, R_current, intrLvl, DM_weights[DM_type], frontier_filter, DT_presentation) 
							# induce a set of rules and sample a unique rule, if a unique solution was not found
							rules = get_rules(gateway)
							# check if at least a rule characterizing 'good' portfolios was induced...
							for rule in rules:
								if 'good' in rule:
									rules_not_induced = False	
							if unique_sol_found:
								print("Unique solution found...")
								rules_not_induced = False	
						if not unique_sol_found:
							selected_rule = sample_rule(rules, DT_presentation) 
						#print(selected_rule)
						# save the simulation inputs
						results_optim["instNum"].append(instNum)
						results_optim["dm_type"].append(DM_type)
						results_optim["intrLvl"].append(intrLvl)
						results_optim["frontier_filter"].append(frontier_filter)
						results_optim["simulation"].append(exc)
						results_optim["dt_type"].append(DT_presentation)
						# save the simulation outputs
						results_optim["current_population"].append([ind.tolist() for ind in R_current])
						results_optim["optimTime"].append(total_CPUTime)
						results_optim["bestAlternative"].append(bestAlternative.tolist())
						fit_best = get_fitness_individuals([bestAlternative], rets_is, I_is)
						fit_best = fit_best[0]
						results_optim["prefFunc"].append(get_pref_func_rebal(np.array(fit_best),DM_weights[DM_type]))
						results_optim["selected_rule"].append(selected_rule)
				# save the final file
				df_optim = pd.DataFrame.from_dict(results_optim)
				df_optim.to_csv(filepath + "results_optim.csv")
				
				# open the file
				simsDB = pd.read_csv(filepath + "results_optim.csv") 
				results_oos = dict([])
				# simulation inputs
				results_oos["t"] = [] # out-of-sample time index
				results_oos["instNum"] = [] # 1,...,6
				results_oos["dm_type"] = [] # TE, ER or EQ
				results_oos["intrLvl"] = [] # interaction level: 1,2,...,intrLvl_max
				results_oos["frontier_filter"] = [] # closer or farther
				results_oos["simulation"] = [] # simulation number: exec
				results_oos["dt_type"] = [] # visual, par_ret, nonpar_ret
				# simulation outputs
				results_oos["cumRet"] = [] 
				results_oos["mean(TE)"] = [] 
				results_oos["max(TE)"] = [] 
				results_oos["mean(ER)"] = [] 
				results_oos["max(-ER)"] = [] 
				results_oos["feasibility"] = [] 
				results_oos["preferenceFunc"] = [] 
				results_oos["portfolio"] = []
				results_oos["portfolio_rebalanced"] = []
				# start the out-of-sample simulation which will be saved in results_oos.csv
				# evaluate the best alternatives obtained in each interaction in the out-of-sample period
				# we only consider portfolios generated after the first interaction: intrLvl > 1
				for DT_presentation in DT_presentations:
					print("DT_presentation (out-of-sample): " + str(DT_presentation))
					for intrLvl in range(2,intrLvl_max+1): # loop over all interaction levels
						current_rebalances_available = max_num_of_rebals
						#print("intrLvl = " + str(intrLvl))
						# get best portfolio composition
						query = (simsDB["intrLvl"] == intrLvl) & (simsDB["dt_type"] == DT_presentation) 
						this_row = simsDB[query]
						portfolio = this_row["bestAlternative"].values[0].replace("[","").replace("]","")
						portfolio = np.fromstring(portfolio,dtype=float,sep = ',')
						# get the population S(intrLvl)
						R = this_row["current_population"].values[0].replace("[","").replace("]","")
						R_current = np.reshape(np.fromstring(R,dtype=float,sep = ','),(nIndividuals,rets.shape[0]))
						# get rule that was obtained from  S(intrLvl-1) 
						query = (simsDB["intrLvl"] == intrLvl-1) & (simsDB["dt_type"] == DT_presentation)
						this_row = simsDB[query]
						rule = json.loads(this_row["selected_rule"].values[0].replace("'", "\""))
						cumRet = 1
						for t in range(T,T+oosT): # loop over all the out-of-sample period
							# simulation inputs
							results_oos["t"].append(t) # t
							results_oos["instNum"].append(instNum)
							results_oos["dm_type"].append(DM_type)
							results_oos["intrLvl"].append(intrLvl)
							results_oos["frontier_filter"].append(frontier_filter)
							results_oos["simulation"].append(exc)
							results_oos["dt_type"].append(DT_presentation)
							t_min = t-T+1
							t_max = t+1
							rets_oos = rets[:,t_min:t_max] # assets returns in-sample
							I_oos = I[t_min:t_max] # assets returns in-sample
							# verify if the current portfolio is still feasible
							fit = get_fitness_individuals([portfolio], rets_oos, I_oos)
							fit = fit[0]
							drsa_eval = get_drsa_eval([portfolio], rets_oos, I_oos, DT_presentation)
							feasibility = get_feasibility_individuals(np.array([portfolio]), rule, drsa_eval)
							# get DM preference in t
							prefFunc_oos= get_pref_func_rebal(np.array(fit), DM_weights[DM_type])
							#print(prefFunc_oos)
							#cumulative return
							cumRet = cumRet*(1+(portfolio @ rets_oos))
							# simulation outputs
							results_oos["cumRet"].append(cumRet)
							results_oos["mean(TE)"].append(fit[0]) 
							results_oos["max(TE)"].append(fit[1])  
							results_oos["mean(ER)"].append(fit[2]) 
							results_oos["max(-ER)"].append(fit[3]) 
							results_oos["feasibility"].append(feasibility[0])
							results_oos["preferenceFunc"].append(prefFunc_oos) 
							results_oos["portfolio"].append(portfolio)
							rebalanced = 0
							# if infeasible and there is still at least one time step available
							if feasibility[0] < 1 and t < T+oosT-1 and current_rebalances_available > 0: # then rebalance the portfolio...
								rebalanced = 1
								#print("reoptimize because infeasible at t = " + str(t) + " (rebalances left: " + str(current_rebalances_available-1) + ")")
								# reoptimize using the rule
								CCEF, R_current, CPUTime = compute_CCEF_rule_guided(nIndividuals, R_current, rule, rets_oos, I_oos, k,  DT_presentation, max_time_rebal)
								portfolio = get_bestAlternative_rebal(rets_oos, I_oos, rule, DT_size, R_current, frontier_filter, DT_presentation) 
								current_rebalances_available -= 1
							results_oos["portfolio_rebalanced"].append(rebalanced)
				# save the final file
				df_oos = pd.DataFrame.from_dict(results_oos)
				df_oos.to_csv(filepath + "results_oos.csv")
'''
# experiments 2 and 3
instances = [1,2,3,4,5,6]
DM_types = ['TE','ER','EQ'] # 3 options: 1- TE (prefers tracking error), 2- ER (prefers excess return), 3 - EQ (Equilibrated investor - equal weights) 
start_exec = 1
frontier_filters = ['closer']
#DT_presentations = ['visual','nonpar_quant'] # for experiment 2
DT_presentations = ['visual'] # for experiment 3
# GA rebalance configurations
max_time_rebal = max_time_S # max time for portfolio rebalance
budget_time_rebal = 300 #(after single period evaluation we will select the best frontier filter and DT_presentation)
max_num_of_rebals = int(budget_time_rebal/max_time_rebal)

for instNum in instances:
	# 1. Read the data 
	thisInstance = Instance(instNum, 'beasley')
	rets = np.array(thisInstance.returns_i) # assets returns 
	I = np.array(thisInstance.returns_Index) # index returns 
	# set the cardinality
	k = card_inst[instNum]
	for DM_type in DM_types:
		# NSGA-II versus IMO-DRSA experiment (comment the code associated with experiment1 (optimizing with pref info) and experiment2 (dynamic IT))
		print("running experiments using NSGA2 without pref info to be compared with the in-sample results of the select DT approaches and the 'closer' filter...")
		# open the file
		filepath = "results/core_experiments/indtrack" + str(instNum) + "/DM_" + DM_type + "/"
		rets_is = rets[:,0:T] # assets returns in-sample
		I_is = I[0:T] # assets returns in-sample
		results_optim = dict([])
		# simulation inputs
		results_optim["instNum"] = [] # 1,...,6
		results_optim["dm_type"] = [] # TE, ER or EQ
		results_optim["intrLvl"] = [] # interaction level: 1,2,...,intrLvl_max
		results_optim["frontier_filter"] = [] # closer or farther
		results_optim["simulation"] = [] # simulation number: exec
		results_optim["dt_type"] = [] # visual, par_ret, nonpar_ret
		# simulation outputs
		results_optim["current_population"] = [] # current population after the optimization (string)
		results_optim["optimTime"] = [] # time until satisfying the proportion of feasible individuals (gamma) in seconds
		results_optim["bestAlternative"] = [] # best alternative if the DM choose an alternative instead of choosing to continue the interaction
		results_optim["prefFunc"] = [] # evaluation of all the elemnts that constituted DT^y (list)
		results_optim["selected_rule"] = [] # rule selected in the current interaction (string)
		for exc in range(1,numRuns+1): 
			# start the simulation
			R_current = []
			selected_rule = []
			for intrLvl in range(1,intrLvl_max+1): # loop over all interaction levels
				#print("generating S" + str(intrLvl-1))
				if intrLvl == 1: # get the initial set of solutions S^0 + data_table1
					CCEF, R_current, CPUTime = compute_CCEF(nIndividuals, rets_is, I_is, k, max_time_S0)
				else: # generate S(intrLvl-1).txt using selected_rule and R_current
					CCEF, R_current, CPUTime = compute_CCEF_not_rule_guided(nIndividuals, R_current, rets_is, I_is, k, max_time_S)
				total_CPUTime = CPUTime
				# generate data_table(intrLvl)
				bestAlternative = get_best_alternative_nsga2(rets_is, I_is, R_current, DM_weights[DM_type]) 
				# save the simulation inputs
				results_optim["instNum"].append(instNum)
				results_optim["dm_type"].append(DM_type)
				results_optim["intrLvl"].append(intrLvl)
				results_optim["frontier_filter"].append('')
				results_optim["simulation"].append(exc)
				results_optim["dt_type"].append('')
				# save the simulation outputs
				results_optim["current_population"].append([ind.tolist() for ind in R_current])
				results_optim["optimTime"].append(total_CPUTime)
				results_optim["bestAlternative"].append(bestAlternative.tolist())
				fit_best = get_fitness_individuals([bestAlternative], rets_is, I_is)
				fit_best = fit_best[0]
				results_optim["prefFunc"].append(get_pref_func_rebal(np.array(fit_best),DM_weights[DM_type]))
				results_optim["selected_rule"].append(selected_rule)
		# save the final file
		df_optim = pd.DataFrame.from_dict(results_optim)
		df_optim.to_csv(filepath + "results_nsga2_nopref.csv")
		
		'''
		for frontier_filter in frontier_filters:
			print("running the dynamic IT experiments using NSGA2 with pref info for the selected DT approaches and the 'closer' filter in the oos period...")
			for exc in range(start_exec,numRuns+1): 
				# open the file
				filepath = "results/core_experiments/indtrack" + str(instNum) + "/DM_" + DM_type + "/" + frontier_filter +"/sim" + str(exc) + "/"
				print("writing: " + filepath)
				simsDB = pd.read_csv(filepath + "results_optim.csv") 
				results_oos = dict([])
				# simulation inputs
				results_oos["t"] = [] # out-of-sample time index
				results_oos["instNum"] = [] # 1,...,6
				results_oos["dm_type"] = [] # TE, ER or EQ
				results_oos["intrLvl"] = [] # interaction level: 1,2,...,intrLvl_max
				results_oos["frontier_filter"] = [] # closer or farther
				results_oos["simulation"] = [] # simulation number: exec
				results_oos["dt_type"] = [] # visual, par_ret, nonpar_ret
				# simulation outputs
				results_oos["cumRet"] = [] 
				results_oos["mean(TE)"] = [] 
				results_oos["max(TE)"] = [] 
				results_oos["mean(ER)"] = [] 
				results_oos["max(-ER)"] = [] 
				results_oos["feasibility"] = [] 
				results_oos["preferenceFunc"] = [] 
				results_oos["portfolio"] = []
				results_oos["portfolio_rebalanced"] = []
				# start the out-of-sample simulation which will be saved in results_oos.csv
				# evaluate the best alternatives obtained in each interaction in the out-of-sample period
				# we only consider portfolios generated after the first interaction: intrLvl > 1
				for DT_presentation in DT_presentations:
					print("DT_presentation (out-of-sample): " + str(DT_presentation))
					for intrLvl in range(2,intrLvl_max+1): # loop over all interaction levels
						current_rebalances_available = max_num_of_rebals
						#print("intrLvl = " + str(intrLvl))
						# get best portfolio composition
						query = (simsDB["intrLvl"] == intrLvl) & (simsDB["dt_type"] == DT_presentation) 
						this_row = simsDB[query]
						portfolio = this_row["bestAlternative"].values[0].replace("[","").replace("]","")
						portfolio = np.fromstring(portfolio,dtype=float,sep = ',')
						# get the population S(intrLvl)
						R = this_row["current_population"].values[0].replace("[","").replace("]","")
						R_current = np.reshape(np.fromstring(R,dtype=float,sep = ','),(nIndividuals,rets.shape[0]))
						# get rule that was obtained from  S(intrLvl-1) 
						query = (simsDB["intrLvl"] == intrLvl-1) & (simsDB["dt_type"] == DT_presentation)
						this_row = simsDB[query]
						rule = json.loads(this_row["selected_rule"].values[0].replace("'", "\""))
						cumRet = 1
						for t in range(T,T+oosT): # loop over all the out-of-sample period
							# simulation inputs
							results_oos["t"].append(t) # t
							results_oos["instNum"].append(instNum)
							results_oos["dm_type"].append(DM_type)
							results_oos["intrLvl"].append(intrLvl)
							results_oos["frontier_filter"].append(frontier_filter)
							results_oos["simulation"].append(exc)
							results_oos["dt_type"].append(DT_presentation)
							t_min = t-T+1
							t_max = t+1
							rets_oos = rets[:,t_min:t_max] # assets returns in-sample
							I_oos = I[t_min:t_max] # assets returns in-sample
							# verify if the current portfolio is still feasible
							fit = get_fitness_individuals([portfolio], rets_oos, I_oos)
							fit = fit[0]
							drsa_eval = get_drsa_eval([portfolio], rets_oos, I_oos, DT_presentation)
							feasibility = get_feasibility_individuals(np.array([portfolio]), rule, drsa_eval)
							# get DM preference in t
							prefFunc_oos= get_pref_func_rebal(np.array(fit), DM_weights[DM_type])
							#print(prefFunc_oos)
							#cumulative return
							cumRet = cumRet*(1+(portfolio @ rets_oos))
							# simulation outputs
							results_oos["cumRet"].append(cumRet)
							results_oos["mean(TE)"].append(fit[0]) 
							results_oos["max(TE)"].append(fit[1])  
							results_oos["mean(ER)"].append(fit[2]) 
							results_oos["max(-ER)"].append(fit[3]) 
							results_oos["feasibility"].append(feasibility[0])
							results_oos["preferenceFunc"].append(prefFunc_oos) 
							results_oos["portfolio"].append(portfolio)
							rebalanced = 0
							# if infeasible and there is still at least one time step available
							if feasibility[0] < 1 and t < T+oosT-1 and current_rebalances_available > 0: # then rebalance the portfolio...
								rebalanced = 1
								#print("reoptimize because infeasible at t = " + str(t) + " (rebalances left: " + str(current_rebalances_available-1) + ")")
								# reoptimize using the rule
								CCEF, R_current, CPUTime = compute_CCEF_rule_guided(nIndividuals, R_current, rule, rets_oos, I_oos, k,  DT_presentation, max_time_rebal)
								portfolio = get_bestAlternative_rebal(rets_oos, I_oos, rule, DT_size, R_current, frontier_filter, DT_presentation) 
								current_rebalances_available -= 1
							results_oos["portfolio_rebalanced"].append(rebalanced)
				# save the final file
				df_oos = pd.DataFrame.from_dict(results_oos)
				df_oos.to_csv(filepath + "results_oos2.csv")
		'''
				














