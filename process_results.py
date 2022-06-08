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
instances = [1,2,3]
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

sns.set(rc={'figure.figsize':(9.7,6.27)})

for instNum in instances:
	# aggregate data for a given instance
	insampleDB = []
	oosampleDB = []
	first_batch = True
	for DM_type in DM_types:
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
	# generate figs
	oosfig_filepath = "results/figures/oosample/" 
	insfig_filepath = "results/figures/insample/"
	if not os.path.exists(oosfig_filepath):
		os.makedirs(oosfig_filepath) # create this path if it not exists
		os.makedirs(insfig_filepath) # create this path if it not exists
	# oos figures
	print("generating figures in " + oosfig_filepath + " associated with indtrack" + str(instNum))
	sns.set_theme(style="ticks", font_scale=1.0)
	#sns_plot = sns.boxplot(data=oosampleDB, x="dm_type", y = "preferenceFunc", hue = "dt_type") col = "intrLvl"
	sns_plot = sns.catplot(data=oosampleDB, x="dt_type", y = "feasibility", hue = "intrLvl", row = "dm_type", col = "frontier_filter", kind = 'bar', sharey = False) 
	#sns_plot.figure.savefig(oosfig_filepath + "indtrack" + str(instNum) + "DMpref.png")
	sns_plot.savefig(oosfig_filepath + "indtrack" + str(instNum) + "DMpref.png")
	#sns_plot.clf()
