import numpy as np
import itertools as it
import cplex as cp

def solve_MV(r_exp, r_var, n_points):

	lmbd = 0 #risk_avrs_coef
	UEF = dict([])
	delta = 1/n_points
	for count in range(n_points):

		print("calculating point " + str(count))

		solution_time = None
		f_val = None
		weights = []
		k = len(r_exp)
		nAssets = k # nAssets + 1 (stock index)
		ub,lb=1,0 

		# compute model parameters
		H=np.zeros([2*nAssets,2*nAssets])
		G=np.zeros([2*nAssets])
	    # calculate matrix H
		for i in range(nAssets):
			for j in range(nAssets):
				H[i,j] = lmbd*r_var[i,j]
	    #calculate matrix gF
		for i in range(nAssets):
			G[i] = (1-lmbd)*r_exp[i]
		
		h=H.tolist()
		g=G.tolist()

		gg=[(i,-g[i])for i in range(nAssets)]+[(i,0) for i in range(nAssets,2*nAssets)]
		square=[(i,j,h[i][j]) for i,j in it.product(range(2*nAssets),range(2*nAssets))]


		c = cp.Cplex()
		var=["w"+ str(i) for i in range(nAssets)]+["z"+ str(i) for i in range(nAssets)]
		c.parameters.timelimit.set(3600)

		indices=c.variables.add(names=var)
		c.variables.set_types((i,c.variables.type.binary) for i in range(nAssets,2*(nAssets)))
		c.objective.set_sense(c.objective.sense.minimize)
		c.objective.set_quadratic_coefficients(square)
		c.objective.set_linear(gg)
		rows=[[var,[1 for i in range(nAssets)]+[0 for i in range(nAssets)]],                  #Porfolio allocation constraint
	          [var,[0 for i in range(nAssets)]+[1 for i in range(nAssets)]]]+[                 #cardinality constraint
	          [var,[0 for i in range(j)]+[1]+[0 for i in range(nAssets-1)]+[-lb]+[0 for i in range(nAssets-1-j)]]for j in range(nAssets)]+[          #lower bound
	          [var,[0 for i in range(j)]+[1]+[0 for i in range(nAssets-1)]+[-ub]+[0 for i in range(nAssets-1-j)]]for j in range(nAssets)]            #upper bound
	 
		indices=c.linear_constraints.add(lin_expr=rows, 
	                                     senses=["E","L"]+["G" for i in range(nAssets)]+["L" for i in range(nAssets)], 
	                                     rhs=[1,k]+[0 for i in range(2*(nAssets))],
	                                     names=["c"+str(i) for i in range(2*(nAssets)+2)])

		c.parameters.mip.display.set(0) # do not display soler outputs
		c.set_log_stream(None)
		c.set_error_stream(None)
		c.set_warning_stream(None)
		c.set_results_stream(None)

		t0=c.get_time()
		c.solve()
		t1=c.get_time()
		solution_time=t1-t0
		f_val=c.solution.get_objective_value()
		x=c.solution.get_values()

		# get mean ret and var
		UEF[lmbd] = [get_mean(x[0:nAssets], r_exp), get_var(x[0:nAssets], r_var)]
		lmbd += delta

	return UEF

def get_mean(weights, r_exp):
	return np.dot(weights, r_exp)

def get_var(weights, r_var):
	first_mult = weights @ r_var
	return first_mult @ weights


