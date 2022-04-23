import os

def write_UEF(instNum, UEF):

	filepath = "results/UEFs/mean_var/" 
	if not os.path.exists(filepath):
		os.makedirs(filepath) # create this path if it not exists

	'''
	write check.txt to notify the metheur_experiements_core.py algorithm that the results
	for the selected_metaheuristic, model, source, instnum and k has already been computed!
	'''
	filehandle = open(filepath + "check.txt", 'w')  

	filepath = filepath + "indtrack" + str(instNum) + ".txt"

	try:
		filehandle = open(filepath, 'w')
		for lmbd in UEF.keys():
			r_exp, r_var = UEF[lmbd]
			filehandle.write(str(r_exp) + '\t' + str(r_var) + '\n')

	except FileNotFoundError:
		print(filepath + " could not be written")

def write_CCEF(instNum, R, exc, intrLvl, DM_type):

	filepath = "samples/indtrack" + str(instNum) + "/" + str(DM_type) + "/exec" + str(exc) + "/S" + str(intrLvl-1) + ".txt"
	# write for risk-averse DM
	print("write in " + filepath)
	filehandle = open(filepath, 'w')
	for ind in range(len(R)):
		portfolio_str = str([r for r in R[ind]]).replace('[','')
		if ind == len(R) - 1:
			portfolio_str = portfolio_str.replace(']','')
		filehandle.write(portfolio_str)

def write_initial_CCEF(instNum, R, exc):

	filepath1 = "samples/indtrack" + str(instNum) + "/risk-averse/exec" + str(exc) + "/"
	filepath2 = "samples/indtrack" + str(instNum) + "/risk-prone/exec" + str(exc) + "/" 
	if not os.path.exists(filepath1):
		os.makedirs(filepath1) # create this path if it not exists
	if not os.path.exists(filepath2):
		os.makedirs(filepath2) # create this path if it not exists

	filepath1 = filepath1 + "S0.txt"
	filepath2 = filepath2 + "S0.txt"

	# write for risk-averse DM
	print("write in " + filepath1)
	filehandle1 = open(filepath1, 'w')
	for ind in range(len(R)):
		portfolio_str = str([r for r in R[ind]]).replace('[','')
		if ind == len(R) - 1:
			portfolio_str = portfolio_str.replace(']','')
		filehandle1.write(portfolio_str)

	# write for risk-prone DM
	print("write in " + filepath2)
	filehandle2 = open(filepath2, 'w')
	for ind in range(len(R)):
		portfolio_str = str([r for r in R[ind]]).replace('[','')
		if ind == len(R) - 1:
			portfolio_str = portfolio_str.replace(']','')
		filehandle2.write(portfolio_str)
