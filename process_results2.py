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
frontier_filters = ['closer'] # best filter
DT_presentations = ['visual','nonpar_quant'] # best DT_presentations: nonpar_quant(non-parametric quantile), visual
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
	insampleDBNSGA2 = []
	first_batch = True
	for DM_type in DM_types:
		for frontier_filter in frontier_filters:
			for exc in range(1,numRuns+1):
				filepath = "results/core_experiments/indtrack" + str(instNum) + "/DM_" + DM_type + "/" + frontier_filter +"/sim" + str(exc) + "/"
				inputDBinsample = pd.read_csv(filepath + "results_optim.csv") 
				inputDBoosample = pd.read_csv(filepath + "results_oos.csv") 
				inputDBoosample2 = pd.read_csv(filepath + "results_oos2.csv") 
				# adjust oosample DBS
				inputDBoosample["Dynamic"] = 'No'
				inputDBoosample2["Dynamic"] = 'Yes'
				if first_batch:
					insampleDB = inputDBinsample
					oosampleDB = inputDBoosample
					oosampleDB2= inputDBoosample2
					first_batch = False
				else:
					insampleDB = pd.concat([insampleDB, inputDBinsample])
					oosampleDB = pd.concat([oosampleDB, inputDBoosample])
					oosampleDB = pd.concat([oosampleDB, inputDBoosample2])
		# add results associated with NSGA2 without preference information
		filepath = "results/core_experiments/indtrack" + str(instNum) + "/DM_" + DM_type + "/"
		inputDBinsample = pd.read_csv(filepath + "results_nsga2_nopref.csv") 
		# apply transformations
		inputDBinsample['dt_type'] = 'none' 
		inputDBinsample['frontier_filter'] = 'none' 
		insampleDB = pd.concat([insampleDB, inputDBinsample])
			
	# apply column name corrections
	oosampleDB = oosampleDB.rename(columns={"preferenceFunc":"prefFunc"})
	oosampleDB2 = oosampleDB2.rename(columns={"preferenceFunc":"prefFunc"})
	# apply transformations
	oosampleDB.loc[oosampleDB['feasibility'] == 1, 'feasibility'] = 0 
	# apply queries
	# adjust oosample DB
	query = (oosampleDB["frontier_filter"] == 'closer') & ((oosampleDB["dt_type"] == 'visual') | (oosampleDB["dt_type"] == 'nonpar_quant'))
	oosampleDB = oosampleDB[query]
	# adjust insample DB
	query = ((insampleDB["frontier_filter"] == 'closer') | (insampleDB["frontier_filter"] == 'none')) & ((insampleDB["dt_type"] == 'visual') | (insampleDB["dt_type"] == 'none') )
	insampleDB = insampleDB[query]
	#query2 = insampleDB["intrLvl"] > 1
	#insampleDB = insampleDB[query2]

	# generate figs
	oosfig_filepath = "results/figures/oosample2/" 
	oosfig_best_filepath = "results/figures/oosample2/best/" 
	insfig_filepath = "results/figures/nsga2/"
	if not os.path.exists(oosfig_filepath):
		os.makedirs(oosfig_filepath) # create this path if it not exists
		os.makedirs(insfig_filepath) # create this path if it not exists
	if not os.path.exists(oosfig_best_filepath):
		os.makedirs(oosfig_best_filepath)

	# set fig configs
	sns.set(rc={'figure.figsize':(9.7,6.27)})
	sns.set_theme(style="ticks", font_scale=1.28)


	# Comparing NSGA2 with and without preference info
	sns.set_theme(style="ticks", font_scale=1.4)
	sns_plot = sns.catplot(data=insampleDB, x="dt_type", y = "prefFunc", hue = "intrLvl", col = "dm_type", kind = 'box', sharey = False,
							order = ['visual', 'none'], col_order = ['TE', 'ER', 'EQ'], height=5.00, aspect=1.0) 
	sns_plot.savefig(insfig_filepath + "indtrack" + str(instNum) + "DMpref_nsga2.png")

	# oosample - comparing single and multiperiod optimization using preference information....
	sns.set_theme(style="ticks", font_scale=1.8)
	sns_plot = sns.relplot(data = oosampleDB, x="t", y="prefFunc", hue="Dynamic", style = 'dt_type', row = 'dm_type', col = 'intrLvl', kind = 'line',
							facet_kws = {'sharey': False, 'sharex': True}, height=5.00, aspect=1.0)
	sns_plot.savefig(oosfig_filepath + "indtrack" + str(instNum) + "DMpref_oos_dynamic.png")

	#comparing the best method in oos dynamic
	query = (oosampleDB["dt_type"] == 'visual') & (oosampleDB["Dynamic"] == 'Yes')
	oosampleDB_best = oosampleDB[query]
	sns.set_theme(style="ticks", font_scale=1.8)
	sns_plot = sns.relplot(data = oosampleDB_best, x="t", y="prefFunc", hue='intrLvl', style = 'dt_type', col = 'dm_type', kind = 'line',
							facet_kws = {'sharey': False, 'sharex': True}, height=5.00, aspect=1.0)
	sns_plot.savefig(oosfig_best_filepath + "indtrack" + str(instNum) + "DMpref_oos_dynamic.png")

	plt.close('all')



#def to_1D(series):
#	return pd.Series([x for _list in series for x in _list])


