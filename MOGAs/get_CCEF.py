import numpy as np
from models.card_mean_var import solve_MV
from MOGAs.nsga2.nsga2 import *
from MOGAs.nsga2.nsga2_with_init_pop import *
from MOGAs.nsga2_rule_guided.nsga2_rule_guided import *
#from result_analysis.plot_UEF import *
import time

def compute_CCEF(nind, rets, I, k, max_time):

	#print("get initial cardinality constrained Eff. Frontier using NSGA-II")
	t_start = time.time() # init time
	CCEF, R = nsga2(nind, rets, I, k, max_time)
	CPUTime = time.time() - t_start # total time

	return CCEF, R, CPUTime

def compute_CCEF_rule_guided(nind, X, rule, rets, I, k, DT_presentation, max_time):

	#print("get cardinality constrained Eff. Frontier using NSGA-II guided by DRSA rules")
	t_start = time.time() # init time
	CCEF, R, feasibility = nsga2_rule_guided(nind, np.array(X), rule, rets, I, k, DT_presentation, max_time)
	CPUTime = time.time() - t_start # total time
	return CCEF, R, CPUTime

def compute_CCEF_not_rule_guided(nind, X, rets, I, k, max_time):

	#print("get cardinality constrained Eff. Frontier using NSGA-II guided by DRSA rules")
	t_start = time.time() # init time
	CCEF, R = nsga2_with_init_pop(nind, np.array(X), rets, I, k, max_time)
	CPUTime = time.time() - t_start # total time
	return CCEF, R, CPUTime

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
