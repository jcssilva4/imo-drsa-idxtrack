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
# stat tests
from scipy.stats import wilcoxon
from scipy.stats import levene

# model and ga parameters
instances = [1,2,3,4,5,6]
#instances = [1,2]
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

# aggregate data for a given instance
insampleDB = []
#oosampleDB = []
first_batch = True
for instNum in instances:
	print("getting data from instance" + str(instNum))
	for DM_type in DM_types:
		for frontier_filter in frontier_filters:
			for exc in range(1,numRuns+1):
				filepath = "results/core_experiments/indtrack" + str(instNum) + "/DM_" + DM_type + "/" + frontier_filter +"/sim" + str(exc) + "/"
				inputDBinsample = pd.read_csv(filepath + "results_optim.csv") 
				#inputDBoosample = pd.read_csv(filepath + "results_oos.csv") 
				if first_batch:
					insampleDB = inputDBinsample
					#oosampleDB = inputDBoosample
					first_batch = False
				else:
					insampleDB = pd.concat([insampleDB, inputDBinsample])
					#oosampleDB = pd.concat([oosampleDB, inputDBoosample])
		# add results associated with NSGA2 without preference information
		filepath = "results/core_experiments/indtrack" + str(instNum) + "/DM_" + DM_type + "/"
		inputDBinsample = pd.read_csv(filepath + "results_nsga2_nopref.csv") 
		# apply transformations
		inputDBinsample['dt_type'] = 'none' 
		inputDBinsample['frontier_filter'] = 'none' 
		insampleDB = pd.concat([insampleDB, inputDBinsample])
					
	# apply column name corrections
	#oosampleDB = oosampleDB.rename(columns={"preferenceFunc":"prefFunc"})
	# apply transformations
	#oosampleDB.loc[oosampleDB['feasibility'] == 1, 'feasibility'] = 0 

# generate dirs
ins_filepath = "results/stat_tests/insample/"
if not os.path.exists(ins_filepath):
	os.makedirs(ins_filepath) # create this path if it not exists

insampleDB = insampleDB[['instNum','intrLvl','frontier_filter','prefFunc','dt_type','dm_type']]


# Wilcoxon tests: closer x farther filter
results_pvalues = dict([])
results_pvalues['dt_type'] = []
results_pvalues['instNum'] = []
results_pvalues['dm_type'] = []
results_pvalues['intrLvl'] = []
results_pvalues['pvalue'] = []
for dt_type in DT_presentations:
	for instance in instances:
		for dm_type in DM_types:
			for intrLvl in range(2,intrLvl_max+1):
				# 'closer' data
				query = (insampleDB["frontier_filter"] == 'closer') & (insampleDB["dt_type"] == dt_type) & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl) 
				filtered_insampleDB = insampleDB[query] 
				values_closer = filtered_insampleDB['prefFunc'].values
				# 'farther' data
				query = (insampleDB["frontier_filter"] == 'farther') & (insampleDB["dt_type"] == dt_type) & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl) 
				filtered_insampleDB = insampleDB[query] 
				values_farther = filtered_insampleDB['prefFunc'].values
				statistic, pvalue = wilcoxon(values_closer,values_farther)
				#print(pvalue)
				results_pvalues['dt_type'].append(dt_type)
				results_pvalues['instNum'].append(instance)
				results_pvalues['dm_type'].append(dm_type)
				results_pvalues['intrLvl'].append(intrLvl)
				results_pvalues['pvalue'].append(pvalue)

# write the results in tables....
df_pvalues = pd.DataFrame.from_dict(results_pvalues)
f = open(ins_filepath + "closer_farther_pvals.txt", 'w')  
for dt_type in DT_presentations:
	#f = open(ins_filepath + "closer_farther_pvals_" + dt_type + ".txt", 'w')  
	f.write('\\begin{table}[h]\n')
	f.write('\\centering\n')
	f.write('\\label{tab:'+ dt_type + '}\n')
	f.write('\\begin{tabular}{c c cccccc}\n')
	f.write('\\toprule\n')
	f.write('\\multirow{2}{*}{dm\\_type} &\n') 
	f.write('\\multirow{2}{*}{intrLvl} &\n') 
	f.write('\\multicolumn{6}{c}{indtrack}    \\\\ \n')
	f.write('\\cmidrule(lr){3-8} \n') 
	f.write('& & 1  & 2  & 3  & 4  & 5  & 6          \\\\ \n')
	f.write('\\midrule \n')
	for dm_type in DM_types:
		f.write("\\multirow{3}{*}{" + dm_type + "} \n")
		for intrLvl in range(2,intrLvl_max+1):
			f.write('& ' + str(intrLvl))
			for instNum in instances:
				query = (df_pvalues["dt_type"] == dt_type) & (df_pvalues["dm_type"] == dm_type) & (df_pvalues["intrLvl"] == intrLvl) & (df_pvalues["instNum"] == instNum)
				filtered_dfpval = df_pvalues[query]
				pval = filtered_dfpval['pvalue'].values[0]
				pval_string = "{:.5f}".format(pval)
				if pval < 0.05:
					f.write('& \\textbf{' + pval_string + '}')
				else:
					f.write('& ' + pval_string)
			f.write('\\\\ \n')
		f.write('\\addlinespace\n')
	f.write('\\bottomrule\n')
	f.write('\\end{tabular} \n')
	f.write('\\caption{Resultant p-values from the comparison between frontier filters for the \'' + dt_type + '\' data table generation method.} \n')
	f.write('\\end{table} \n\n')

	

#wilcoxon tests - difference along interactions...
results_pvalues = dict([])
results_pvalues['dt_type'] = []
results_pvalues['instNum'] = []
results_pvalues['dm_type'] = []
results_pvalues['intrLvl'] = []
results_pvalues['pvalue'] = []
for dt_type in DT_presentations:
	for instance in instances:
		for dm_type in DM_types:
			for intrLvl in range(2,intrLvl_max+1):
				# 'intrLvl-1 data
				query = (insampleDB["frontier_filter"] == 'closer') & (insampleDB["dt_type"] == dt_type) & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl-1) 
				filtered_insampleDB = insampleDB[query] 
				values_pastIntrlvl = filtered_insampleDB['prefFunc'].values
				# intrLvl data
				query = (insampleDB["frontier_filter"] == 'closer') & (insampleDB["dt_type"] == dt_type) & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl) 
				filtered_insampleDB = insampleDB[query] 
				values_currentIntrlvl = filtered_insampleDB['prefFunc'].values
				statistic, pvalue = wilcoxon(values_currentIntrlvl,values_pastIntrlvl)
				results_pvalues['dt_type'].append(dt_type)
				results_pvalues['instNum'].append(instance)
				results_pvalues['dm_type'].append(dm_type)
				results_pvalues['intrLvl'].append(intrLvl)
				results_pvalues['pvalue'].append(pvalue)

# write the results in tables....
df_pvalues = pd.DataFrame.from_dict(results_pvalues)
f = open(ins_filepath + "intrLvl_pvals.txt", 'w')  
for dt_type in DT_presentations:
	#f = open(ins_filepath + "closer_farther_pvals_" + dt_type + ".txt", 'w')  
	f.write('\\begin{table}[h]\n')
	f.write('\\centering\n')
	f.write('\\label{tab:'+ dt_type + '}\n')
	f.write('\\begin{tabular}{c c cccccc}\n')
	f.write('\\toprule\n')
	f.write('\\multirow{2}{*}{dm\\_type} &\n') 
	f.write('\\multirow{2}{*}{(intrLvl-1, intrLvl)} &\n') 
	f.write('\\multicolumn{6}{c}{indtrack}    \\\\ \n')
	f.write('\\cmidrule(lr){3-8} \n') 
	f.write('& & 1  & 2  & 3  & 4  & 5  & 6          \\\\ \n')
	f.write('\\midrule \n')
	for dm_type in DM_types:
		f.write("\\multirow{3}{*}{" + dm_type + "} \n")
		for intrLvl in range(2,intrLvl_max+1):
			f.write('& (' + str(intrLvl-1) + ", " + str(intrLvl) + ")")
			for instNum in instances:
				query = (df_pvalues["dt_type"] == dt_type) & (df_pvalues["dm_type"] == dm_type) & (df_pvalues["intrLvl"] == intrLvl) & (df_pvalues["instNum"] == instNum)
				filtered_dfpval = df_pvalues[query]
				pval = filtered_dfpval['pvalue'].values[0]
				pval_string = "{:.5f}".format(pval)
				if pval < 0.05:
					f.write('& \\textbf{' + pval_string + '}')
				else:
					f.write('& ' + pval_string)
			f.write('\\\\ \n')
		f.write('\\addlinespace\n')
	f.write('\\bottomrule\n')
	f.write('\\end{tabular} \n')
	f.write('\\caption{Resultant p-values from the comparison between interaction levels for the \'' + dt_type + '\' data table generation method.} \n')
	f.write('\\end{table} \n\n')

# NSGA2 x IMO-DRSA (closer filter + visual approach)
# adjust insample DB
query = ((insampleDB["frontier_filter"] == 'closer') | (insampleDB["frontier_filter"] == 'none')) & ((insampleDB["dt_type"] == 'visual') | (insampleDB["dt_type"] == 'none') )
DT_presentations_nsga2 = ['visual','none']
#wilcoxon tests - difference along interactions...
results_pvalues = dict([])
results_pvalues['instNum'] = []
results_pvalues['dm_type'] = []
results_pvalues['intrLvl'] = []
results_pvalues['pvalue'] = []
for instance in instances:
	for dm_type in DM_types:
		for intrLvl in range(2,intrLvl_max+1):
			# imodrsa data
			query = (insampleDB["frontier_filter"] == 'closer') & (insampleDB["dt_type"] == 'visual') & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl) 
			filtered_insampleDB = insampleDB[query] 
			values_visual = filtered_insampleDB['prefFunc'].values
			# nsga2 data
			query = (insampleDB["frontier_filter"] == 'none') & (insampleDB["dt_type"] == 'none') & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl) 
			filtered_insampleDB = insampleDB[query] 
			values_none = filtered_insampleDB['prefFunc'].values
			statistic, pvalue = wilcoxon(values_visual,values_none)
			results_pvalues['instNum'].append(instance)
			results_pvalues['dm_type'].append(dm_type)
			results_pvalues['intrLvl'].append(intrLvl)
			results_pvalues['pvalue'].append(pvalue)

# write the results in tables....
df_pvalues = pd.DataFrame.from_dict(results_pvalues)
f = open(ins_filepath + "nsga2_pvals.txt", 'w')  
#f = open(ins_filepath + "closer_farther_pvals_" + dt_type + ".txt", 'w')  
f.write('\\begin{table}[h]\n')
f.write('\\centering\n')
f.write('\\label{tab:'+ dt_type + '}\n')
f.write('\\begin{tabular}{c c cccccc}\n')
f.write('\\toprule\n')
f.write('\\multirow{2}{*}{dm\\_type} &\n') 
f.write('\\multirow{2}{*}{(intrLvl)} &\n') 
f.write('\\multicolumn{6}{c}{indtrack}    \\\\ \n')
f.write('\\cmidrule(lr){3-8} \n') 
f.write('& & 1  & 2  & 3  & 4  & 5  & 6          \\\\ \n')
f.write('\\midrule \n')
for dm_type in DM_types:
	f.write("\\multirow{3}{*}{" + dm_type + "} \n")
	for intrLvl in range(2,intrLvl_max+1):
		f.write('& '+ str(intrLvl))
		for instNum in instances:
			query = (df_pvalues["dm_type"] == dm_type) & (df_pvalues["intrLvl"] == intrLvl) & (df_pvalues["instNum"] == instNum)
			filtered_dfpval = df_pvalues[query]
			pval = filtered_dfpval['pvalue'].values[0]
			pval_string = "{:.5f}".format(pval)
			if pval < 0.05:
				f.write('& \\textbf{' + pval_string + '}')
			else:
				f.write('& ' + pval_string)
		f.write('\\\\ \n')
	f.write('\\addlinespace\n')
f.write('\\bottomrule\n')
f.write('\\end{tabular} \n')
f.write('\\caption{Resultant p-values from the comparison between NSGA-II and IMO-DRSA with the closer frontier filter and visual approach in each interaction level.} \n')
f.write('\\end{table} \n\n')



# final tests for the best apprach: 'visual'

# Levene tests ('median'): closer x farther filter with 'visual' approach
results_pvalues = dict([])
results_pvalues['instNum'] = []
results_pvalues['dm_type'] = []
results_pvalues['intrLvl'] = []
results_pvalues['pvalue'] = []
results_pvalues['sigma2diff'] = []
dt_type = 'visual'
for instance in instances:
	for dm_type in DM_types:
		for intrLvl in range(2,intrLvl_max+1):
			# 'closer' data
			query = (insampleDB["frontier_filter"] == 'closer') & (insampleDB["dt_type"] == dt_type) & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl) 
			filtered_insampleDB = insampleDB[query] 
			values_closer = filtered_insampleDB['prefFunc'].values
			# 'farther' data
			query = (insampleDB["frontier_filter"] == 'farther') & (insampleDB["dt_type"] == dt_type) & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl) 
			filtered_insampleDB = insampleDB[query] 
			values_farther = filtered_insampleDB['prefFunc'].values
			statistic, pvalue = levene(values_closer,values_farther)
			sigma2diff = np.var(values_closer) - np.var(values_farther)
			#print(pvalue)
			results_pvalues['instNum'].append(instance)
			results_pvalues['dm_type'].append(dm_type)
			results_pvalues['intrLvl'].append(intrLvl)
			results_pvalues['pvalue'].append(pvalue)
			results_pvalues['sigma2diff'].append(sigma2diff)

# write the results in tables....
df_pvalues = pd.DataFrame.from_dict(results_pvalues)
f = open(ins_filepath + "levene_visual_closer_farther_pvals.txt", 'w')  
#f = open(ins_filepath + "closer_farther_pvals_" + dt_type + ".txt", 'w')  
f.write('\\begin{table}[h]\n')
f.write('\\centering\n')
f.write('\\label{tab:'+ dt_type + '}\n')
f.write('\\begin{tabular}{c c cccccc}\n')
f.write('\\toprule\n')
f.write('\\multirow{2}{*}{dm\\_type} &\n') 
f.write('\\multirow{2}{*}{intrLvl} &\n') 
f.write('\\multicolumn{6}{c}{indtrack}    \\\\ \n')
f.write('\\cmidrule(lr){3-8} \n') 
f.write('& & 1  & 2  & 3  & 4  & 5  & 6          \\\\ \n')
f.write('\\midrule \n')
for dm_type in DM_types:
	f.write("\\multirow{3}{*}{" + dm_type + "} \n")
	for intrLvl in range(2,intrLvl_max+1):
		f.write('& ' + str(intrLvl))
		for instNum in instances:
			query = (df_pvalues["dm_type"] == dm_type) & (df_pvalues["intrLvl"] == intrLvl) & (df_pvalues["instNum"] == instNum)
			filtered_dfpval = df_pvalues[query]
			pval = filtered_dfpval['pvalue'].values[0]
			sigma2diff = filtered_dfpval['sigma2diff'].values[0]
			pval_string = "{:.5f}".format(pval)
			sigma2diff_string = "{:.2e}".format(sigma2diff)
			f.write('& ' + sigma2diff_string + ' ')
			if pval < 0.05:
				f.write('(\\textbf{' + pval_string + '})')
			else:
				f.write('(' + pval_string + ')')
		f.write('\\\\ \n')
	f.write('\\addlinespace\n')
f.write('\\bottomrule\n')
f.write('\\end{tabular} \n')
f.write('\\caption{Resultant variance differences and p-values of the Levene test from the comparison between frontier filters for the \'' + dt_type + '\' data table generation method.} \n')
f.write('\\end{table} \n\n')


# NSGA2 x IMO-DRSA (closer filter + visual approach)
# adjust insample DB
query = ((insampleDB["frontier_filter"] == 'closer') | (insampleDB["frontier_filter"] == 'none')) & ((insampleDB["dt_type"] == 'visual') | (insampleDB["dt_type"] == 'none') )
DT_presentations_nsga2 = ['visual','none']
#wilcoxon tests - difference along interactions...
results_pvalues = dict([])
results_pvalues['instNum'] = []
results_pvalues['dm_type'] = []
results_pvalues['intrLvl'] = []
results_pvalues['pvalue'] = []
results_pvalues['sigma2diff'] = []
for instance in instances:
	for dm_type in DM_types:
		for intrLvl in range(2,intrLvl_max+1):
			# imodrsa data
			query = (insampleDB["frontier_filter"] == 'closer') & (insampleDB["dt_type"] == 'visual') & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl) 
			filtered_insampleDB = insampleDB[query] 
			values_visual = filtered_insampleDB['prefFunc'].values
			# nsga2 data
			query = (insampleDB["frontier_filter"] == 'none') & (insampleDB["dt_type"] == 'none') & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl) 
			filtered_insampleDB = insampleDB[query] 
			values_none = filtered_insampleDB['prefFunc'].values
			statistic, pvalue = levene(values_visual,values_none)
			sigma2diff = np.var(values_visual) - np.var(values_none)
			results_pvalues['instNum'].append(instance)
			results_pvalues['dm_type'].append(dm_type)
			results_pvalues['intrLvl'].append(intrLvl)
			results_pvalues['pvalue'].append(pvalue)
			results_pvalues['sigma2diff'].append(sigma2diff)

# write the results in tables....
df_pvalues = pd.DataFrame.from_dict(results_pvalues)
f = open(ins_filepath + "levene_nsga2_pvals.txt", 'w')  
#f = open(ins_filepath + "closer_farther_pvals_" + dt_type + ".txt", 'w')  
f.write('\\begin{table}[h]\n')
f.write('\\centering\n')
f.write('\\label{tab:'+ dt_type + '}\n')
f.write('\\begin{tabular}{c c cccccc}\n')
f.write('\\toprule\n')
f.write('\\multirow{2}{*}{dm\\_type} &\n') 
f.write('\\multirow{2}{*}{intrLvl} &\n') 
f.write('\\multicolumn{6}{c}{indtrack}    \\\\ \n')
f.write('\\cmidrule(lr){3-8} \n') 
f.write('& & 1  & 2  & 3  & 4  & 5  & 6          \\\\ \n')
f.write('\\midrule \n')
for dm_type in DM_types:
	f.write("\\multirow{3}{*}{" + dm_type + "} \n")
	for intrLvl in range(2,intrLvl_max+1):
		f.write('& ' + str(intrLvl))
		for instNum in instances:
			query = (df_pvalues["dm_type"] == dm_type) & (df_pvalues["intrLvl"] == intrLvl) & (df_pvalues["instNum"] == instNum)
			filtered_dfpval = df_pvalues[query]
			pval = filtered_dfpval['pvalue'].values[0]
			sigma2diff = filtered_dfpval['sigma2diff'].values[0]
			pval_string = "{:.5f}".format(pval)
			sigma2diff_string = "{:.2e}".format(sigma2diff)
			f.write('& ' + sigma2diff_string + ' ')
			if pval < 0.05:
				f.write('(\\textbf{' + pval_string + '})')
			else:
				f.write('(' + pval_string + ')')
		f.write('\\\\ \n')
	f.write('\\addlinespace\n')
f.write('\\bottomrule\n')
f.write('\\end{tabular} \n')
f.write('\\caption{Resultant variance differences and p-values of the Levene test from the comparison between NSGA-II and IMO-DRSA combining the \'' + dt_type + '\' data table generation method and the \'closer\' frontier filter.} \n')
f.write('\\end{table} \n\n')

#levene - difference along interactions...
results_pvalues = dict([])
results_pvalues['dt_type'] = []
results_pvalues['instNum'] = []
results_pvalues['dm_type'] = []
results_pvalues['intrLvl'] = []
results_pvalues['pvalue'] = []
results_pvalues['sigma2diff'] = []
dt_type = 'visual'
for instance in instances:
	for dm_type in DM_types:
		for intrLvl in range(2,intrLvl_max+1):
			# 'intrLvl-1 data
			query = (insampleDB["frontier_filter"] == 'closer') & (insampleDB["dt_type"] == dt_type) & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl-1) 
			filtered_insampleDB = insampleDB[query] 
			values_pastIntrlvl = filtered_insampleDB['prefFunc'].values
			# intrLvl data
			query = (insampleDB["frontier_filter"] == 'closer') & (insampleDB["dt_type"] == dt_type) & (insampleDB["dm_type"] == dm_type) & (insampleDB["instNum"] == instance) & (insampleDB["intrLvl"] == intrLvl) 
			filtered_insampleDB = insampleDB[query] 
			values_currentIntrlvl = filtered_insampleDB['prefFunc'].values
			statistic, pvalue = levene(values_currentIntrlvl,values_pastIntrlvl)
			sigma2diff = np.var(values_currentIntrlvl) - np.var(values_pastIntrlvl)
			results_pvalues['dt_type'].append(dt_type)
			results_pvalues['instNum'].append(instance)
			results_pvalues['dm_type'].append(dm_type)
			results_pvalues['intrLvl'].append(intrLvl)
			results_pvalues['pvalue'].append(pvalue)
			results_pvalues['sigma2diff'].append(sigma2diff)

# write the results in tables....
df_pvalues = pd.DataFrame.from_dict(results_pvalues)
f = open(ins_filepath + "levene_intrLvl_pvals.txt", 'w')  
#f = open(ins_filepath + "closer_farther_pvals_" + dt_type + ".txt", 'w')  
f.write('\\begin{table}[h]\n')
f.write('\\centering\n')
f.write('\\label{tab:'+ dt_type + '}\n')
f.write('\\begin{tabular}{c c cccccc}\n')
f.write('\\toprule\n')
f.write('\\multirow{2}{*}{dm\\_type} &\n') 
f.write('\\multirow{2}{*}{(intrLvl-1, intrLvl)} &\n') 
f.write('\\multicolumn{6}{c}{indtrack}    \\\\ \n')
f.write('\\cmidrule(lr){3-8} \n') 
f.write('& & 1  & 2  & 3  & 4  & 5  & 6          \\\\ \n')
f.write('\\midrule \n')
for dm_type in DM_types:
	f.write("\\multirow{3}{*}{" + dm_type + "} \n")
	for intrLvl in range(2,intrLvl_max+1):
		f.write('& (' + str(intrLvl-1) + ", " + str(intrLvl) + ")")
		for instNum in instances:
			query = (df_pvalues["dt_type"] == dt_type) & (df_pvalues["dm_type"] == dm_type) & (df_pvalues["intrLvl"] == intrLvl) & (df_pvalues["instNum"] == instNum)
			filtered_dfpval = df_pvalues[query]
			pval = filtered_dfpval['pvalue'].values[0]
			sigma2diff = filtered_dfpval['sigma2diff'].values[0]
			pval_string = "{:.5f}".format(pval)
			sigma2diff_string = "{:.2e}".format(sigma2diff)
			f.write('& ' + sigma2diff_string + ' ')
			if pval < 0.05:
				f.write('(\\textbf{' + pval_string + '})')
			else:
				f.write('(' + pval_string + ')')
		f.write('\\\\ \n')
	f.write('\\addlinespace\n')
f.write('\\bottomrule\n')
f.write('\\end{tabular} \n')
f.write('\\caption{Resultant variance differences and p-values of the Levene test from the comparison between interaction levels for the \'' + dt_type + '\' data table generation method.} \n')
f.write('\\end{table} \n\n')