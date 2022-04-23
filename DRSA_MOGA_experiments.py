from data.DBInstances import Instance
from results.write_results import *
#from get_solution import get_solution
from MOGAs.get_portfolio import get_solution_approx
#from result_analysis.read_cplex_result_file import get_f_val
from MOGAs.get_CCEF import compute_CCEF

selected_instances = dict([])
selected_instances['beasley'] = [5] # change this parameter for each cloud node
k = 10 # select cardinalities
#T = 145 # in-sample period
T = 290
num_executions = 30 # number of runs per heuristic


for instNum in selected_instances['beasley']: # loop over all instances of this source
	
	# 1. Read the data 
	thisInstance = Instance(instNum, 'beasley')

	'''
	filepath = "results/nsga2/" + "indtrack" + str(instNum) + "/check.txt"

	try:
		# check if there are any rules
		print("\n\nchecking: " + filepath)
		#filehandle = open(filepath, 'r')
		#print("results were already computed for indtrack" + str(instNum) + "\n")

		# if there is any rule file, then read the last nsga2 solution file and initialize nsga2 from this solution

	except FileNotFoundError: # if there are no rules, start nsga2 from a random population
	'''

	print("set the optimization problem for indtrack" + str(instNum) + "\n")
	print("model: cardinality constrained bi-objective mean-variance")


	# UEF 
	CCEF = compute_CCEF(thisInstance.returns_i, k, T, instNum)


	# read indtrack_x UEF and plot it with the current nsga2 solution
	# compare_CCEF_UEF(instNum, CCEF)

	# Write frontier
	#write_CCEF(instNum, CCEF)



