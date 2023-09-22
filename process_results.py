from data.DBInstances import Instance
from data_tables import *
# from get_solution import get_solution
from MOGAs.get_portfolio import get_solution_approx
# from result_analysis.read_cplex_result_file import get_f_val
from MOGAs.get_CCEF import *
from rulegen_test import *
import json
# plots
import matplotlib.pyplot as plt 
import seaborn as sns

# model and ga parameters
instances = [1,2,3,4,5,6]
#instances = [1]
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

# simulated IMO-DRSA paramaters
intrLvl_max = 4 # maximum interaction level (lvl1: init sol, lvl2: sol after first interaction,....)
DM_types = ['TE','ER','EQ'] # 3 options: 1- TE (prefers tracking error), 2- ER (preferes excess return), 3 - EQ (Equilibrated investor - equal weights) 
DM_weights = dict([])
DM_weights['TE'] = [0.6, 0.3, 0.05, 0.05]
DM_weights['ER'] = [0.05, 0.05, 0.6, 0.3]
DM_weights['EQ'] = [0.25, 0.25, 0.25, 0.25]
frontier_filters = ['closer','farther']
DT_presentations = ['visual','par_quant','nonpar_quant'] # par_quant (parametric quantile), nonpar_quant(non-parametric quantile), visual
numRuns = 30 # number of simulations per each DM type, instance, and frontier filter (each simulation consists in an interaction process from 1 to intrLvl_max).
DT_size = 3 # Num. of decision examples per data table = 2*DT_size
# GA optimization time
max_time_S0 = 5 # max_time for nsga2
max_time_S = 10 # max_time for rule-guided nsga2
# GA rebalance configurations
max_time_rebal = 2 # max time for portfolio rebalance
budget_time_rebal = 30
max_num_of_rebals = int(budget_time_rebal/max_time_rebal)
# GA population: num of individuals 
nIndividuals = 100

for instNum in instances:
	print("generating plots for instance: " + str(instNum))
	# aggregate data for a given instance
	insampleDB = []
	oosampleDB = []
	first_batch = True
	for DM_type in DM_types:
		print("DM_type: " + DM_type)
		for frontier_filter in frontier_filters:
			for exc in range(1,numRuns+1):
				filepath = "results/core_experiments/indtrack" + str(instNum) + "/DM_" + DM_type + "/" + frontier_filter +"/sim" + str(exc) + "/"
				inputDBinsample = pd.read_csv(filepath + "results_optim.csv") 
				inputDBoosample = pd.read_csv(filepath + "results_oos.csv") 
				if first_batch:
					insampleDB = inputDBinsample
					oosampleDB = inputDBoosample
					first_batch = False
				else:
					insampleDB = pd.concat([insampleDB, inputDBinsample])
					oosampleDB = pd.concat([oosampleDB, inputDBoosample])
					
	# apply column name corrections
	oosampleDB = oosampleDB.rename(columns={"preferenceFunc":"prefFunc"})
	# apply transformations
	oosampleDB.loc[oosampleDB['feasibility'] == 1, 'feasibility'] = 0 

	# generate figs
	oosfig_filepath = "results/figures/oosample/" 
	insfig_filepath = "results/figures/insample/"
	oosfig_best_filepath = oosfig_filepath+'best_front_filter/'
	insfig_best_filepath = insfig_filepath+'best_front_filter/'
	if not os.path.exists(oosfig_filepath) or not os.path.exists(oosfig_best_filepath):
		os.makedirs(oosfig_filepath) # create this path if it not exists
		os.makedirs(insfig_filepath) # create this path if it not exists
		os.makedirs(oosfig_best_filepath) # create this path if it not exists
		os.makedirs(insfig_best_filepath) # create this path if it not exists

	# set fig configs
	#sns.set(rc={'figure.figsize':(20.7,9.27)})
	sns.set_theme(style="ticks", font_scale=1.4)

	# insample figures
	print('insample fig...')
	sns_plot = sns.catplot(data=insampleDB, x="dt_type", y = "prefFunc", hue = "intrLvl", row = "frontier_filter", col = "dm_type", kind = 'box', sharey = False,
							order = ['visual','par_quant','nonpar_quant'], col_order = ['TE', 'ER', 'EQ'], height=5.00, aspect=1.0) 
	sns_plot.savefig(insfig_filepath + "indtrack" + str(instNum) + "DMpref.png")
	# oosample figures
	sns.set_theme(style="ticks", font_scale=1.8)
	print('oosample fig...')
	sns_plot = sns.relplot(data = oosampleDB, x="t", y="prefFunc", hue="frontier_filter", style = 'dt_type', row = 'dm_type', col = 'intrLvl', kind = 'line',
							facet_kws = {'sharey': False, 'sharex': True}, height=5.00, aspect=1.0)
	sns_plot.savefig(oosfig_filepath + "indtrack" + str(instNum) + "DMpref_oos.png")
	

	# the best frontier filter was the 'closer' method, then we continue the evaluations with it...
	print('plot figs only for the closer method...')
	# plot the insample figures, just for the 'closer' method
	# insample figures
	filtered_insampleDB_closer = insampleDB[insampleDB['frontier_filter'] == 'closer']
	#filtered_insampleDB_closer = filtered_insampleDB_closer[filtered_insampleDB_closer['intrLvl']>1]
	sns.set_theme(style="ticks", font_scale=1.4)
	sns_plot = sns.catplot(data=filtered_insampleDB_closer, x="dt_type", y = "prefFunc", hue = "intrLvl",col = "dm_type", kind = 'box', sharey = False,
							order = ['visual','par_quant','nonpar_quant'], col_order = ['TE', 'ER', 'EQ'], height=5.00, aspect=1.0) 
	sns_plot.savefig(insfig_best_filepath + "indtrack" + str(instNum) + "DMpref.png")
	# oosample figures (the best dt_type methods until now: visual and nonpar_quant)
	filtered_oosampleDB_closer = oosampleDB[oosampleDB['frontier_filter'] == 'closer'] 
	#best_dt_query = (filtered_oosampleDB_closer['dt_type'] == 'nonpar_quant') | (filtered_oosampleDB_closer['dt_type'] == 'visual')
	best_dt_query = filtered_oosampleDB_closer['dt_type'] == 'visual'
	filtered_oosampleDB_closer_dt = filtered_oosampleDB_closer[best_dt_query]
	sns.set_theme(style="ticks", font_scale=1.9)
	sns_plot = sns.relplot(data = filtered_oosampleDB_closer_dt, x="t", y="prefFunc", hue='intrLvl', col = 'dm_type', kind = 'line',
							facet_kws = {'sharey': False, 'sharex': True}, height=5.00, aspect=1.0)
	sns_plot.savefig(oosfig_best_filepath + "indtrack" + str(instNum) + "DMpref_oos.png")

	# generating the criteria importance figure
	print('criteria importance...')
	row = 0
	col = 0
	criteria_data = dict([])
	criteria_data["dt_type"] = []
	criteria_data["dm_type"] = []
	criteria_data["intrLvl"] = []
	criteria_data["condition"] = []
	for dt_type in DT_presentations:
		#print("dt_type: " + str(dt_type))
		criteria, trsh1, trsh2 = get_head_DT_presentation(dt_type)
		for dm_type in DM_types:
			for intrLvl in range(2,intrLvl_max+1):
				# query
				query = (insampleDB["frontier_filter"] == 'closer') & (insampleDB["dt_type"] == dt_type) & (insampleDB["dm_type"] == dm_type) & (insampleDB["intrLvl"] == intrLvl-1) 
				filtered_insampleDB = insampleDB[query] 
				# count the num of criteria 
				criteria_dict = {}
				for rule_raw in filtered_insampleDB["selected_rule"]:
					# get the criteria associated with this rule
					rule = json.loads(rule_raw.replace("'", "\""))
					#print(rule)
					count_cond = 0
					cond_dir = rule['cond_dir']
					for j in rule['obj']: # loop over all the associated criteria
						criteria_data["dt_type"].append(dt_type)
						criteria_data["dm_type"].append(dm_type)	
						criteria_data["intrLvl"].append(intrLvl-1)
						criteria_data["condition"].append(criteria[j] + cond_dir[count_cond])
						count_cond += 1
	criteria_df = pd.DataFrame.from_dict(criteria_data)	
	# set fig configs
	#sns.set_theme(style="ticks", font_scale=1.25)
	# generate and save...
	'''
	sns_plot = sns.catplot(data=criteria_df, x="dt_type", hue = "condition", row = "dm_type", col = "intrLvl", kind = 'count', sharey = False, height = 5.7, aspect = 1.0,
								order = ['visual','par_quant','nonpar_quant'], 
								hue_order = ['meanTE>=','meanTE<=','maxTE>=', 'maxTE<=', 'meanER>=', 'meanER<=', 'maxNegER>=', 'maxNegER<=',
								'TE1>=', 'TE1<=', 'TE50>=', 'TE50<=', 'ER1>=', 'ER1<=', 'ER50>=', 'ER50<=', 'ER99>=', 'ER99<=']) 
	'''
	sns.set_theme(style="ticks", font_scale=3.1)
	sns_plot = sns.catplot(data=criteria_df, x="condition", hue = "intrLvl", row = "dm_type", col = "dt_type", kind = 'count', sharex = False, sharey = False,
								col_order = ['visual','par_quant','nonpar_quant'], 
								hue_order = [1,2,3], height=8.00, aspect=1.7)

	
	sns_plot.set_xticklabels(rotation=10)
	sns_plot.savefig(insfig_filepath + "indtrack" + str(instNum) + "obj_importance.png")
	plt.close('all')



#def to_1D(series):
#	return pd.Series([x for _list in series for x in _list])


