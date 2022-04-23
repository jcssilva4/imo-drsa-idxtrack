
''' 
General description of the experiments

Data: {'Beasley': {indtrack1, indtrack2, indtrack3, indtrack4, indtrack5, indtrack6},
	 (Try to get new IBOV data from 2018-2020)}

Models: {Torrubiano2009, Wang2012}

Cardinality constraint: {.05*nAssets, .10*nAssets, .15*nAssets, .20*nAssets}

Results : {Time required to solve (max - 1h), f_val, weight_vector} 

# Total time #
# considering 30 runs performed by the metaheuristic
worst case time per instance per model per cardinality: 2.5h (30*5min = 150min )
worst case time per instance per model: 10h
worst case time per instance: ~80h, which is ~3 days 
using n = nInstances parallel nodes, the worst time required to run all the experiments is = worst case time per instance
use one instance for each metaheuristic

# THE RESULTS FILE #:
directory: /results/cplex/instance/model/
file name: results_k.txt
	first row: 
	Time required to solve (in s)
	second row: 
	f_val
	from the third row onwards : 
	w1
	w2
	.
	.
	.
	w_n
'''

'''
Setting the required google cloud platform machines
# 1. create a VM instance  
# 2. Connect to it clicking in the SSH button at connect column
# 3. sudo apt-get update
# 4. install pip: sudo apt install python3-pip (run those as python3 and pip3, respectively)
# 5. install other required packages associated with your code
# 6. git clone https://github.com/jcssilva4/GRASP-tracking-portfolio.git
# 7. clone instances: click the VM instance and the create a new image from it
# 8. Then go to images and click on the three dots in the final column of the image to create a new instance from it
# 9. after creating all the 7 new VM instances, start the experiments (python 3 cplex_experiments_core.py 'source_name' 'instance_num')

how to run a specific instance.: 

console:

'python 3 cplex_experiments_core.py heuristic_name' # run the script

'ctrl + A + D' # detach from the screen

exit the session

to reopen the sreen use screen -r cplex (or just screen -r)


'''


from data.DBInstances import Instance
from results.write_results import *
#from get_solution import get_solution
from MOGAs.get_portfolio import get_solution_approx
#from result_analysis.read_cplex_result_file import get_f_val
from MOGAs.get_UEF import compute_UEF

selected_instances = dict([])
selected_instances['beasley'] = [5,6] # change this parameter for each cloud node
#selected_instances['beasley'] = [6] # change this parameter for each cloud node
k = 10 # select cardinalities
#T = 145 # in-sample period
oosT = 50 # out-of-sample period
T = 290 - oosT
num_executions = 30 # number of runs per heuristic


for instNum in selected_instances['beasley']: # loop over all instances of this source
	
	# 1. Read the data 
	thisInstance = Instance(instNum, 'beasley')

	filepath = "results/" + "MO_indtrack" + str(instNum) + "/check.txt"

	try:
		# check if the results were already computed for this model, metaheur, source, instnum, k
		print("\n\nchecking: " + filepath)
		filehandle = open(filepath, 'r')
		print("results were already computed for indtrack" + str(instNum) + "\n")

	except FileNotFoundError:
		print("set the optimization problem for _indtrack" + str(instNum) + "\n")
		print("model: cardinality constrained bi-objective mean-variance")

		# 2. Select the model and solve it 
		'''
		res = dict([])
		res['time'], res['f_val'], res ['it'], res['time_to_target'],  res['w'] = get_solution_approx( 
			thisInstance.returns_i, k, T, num_executions)
		'''

		# UEF 
		UEF = compute_UEF(thisInstance.returns_i, k, T)
		# Write UEF
		write_UEF(instNum, UEF)
		'''
		# 3. Write the results in a file
		if not selected_metaheuristic in ['GQ','GQ-MiLS']:
			print("writing GA results")
			write_approx_results_for_k(selected_metaheuristic, source, instNum, model, k, res, num_executions)
		else:
			print("writing GQ results")
			write_GQ_approx_results_for_k(selected_metaheuristic, source, instNum, model, k, res, num_executions)
		'''


