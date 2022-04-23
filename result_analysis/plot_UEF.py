import seaborn as sns
import matplotlib.pyplot as plt


def get_UEF_data(instNum):

	filepath =  "results/UEFs/mean_var/" + "indtrack" + str(instNum) + ".txt"
	r_exp = [] # initialize expected return vec
	r_var = [] # initialize variance vec

	try:

		filelines = []
		filehandle = open(filepath, 'r')
		filelines.extend(filehandle.readlines())
		
		for line in filelines:

			point = line.split('\t')
			r_exp.append(float(point[0]))
			r_var.append(float(point[1].replace('\n', '')))

	except FileNotFoundError:
		print(filepath + " could not be read")

	return r_exp, r_var

def plot_UEfrontier(r_exp, r_var):
	plt.scatter(x = r_var, y = r_exp, s = 6, c = "black")
	plt.xlabel("Variance", fontsize = 18)
	plt.ylabel("Mean", fontsize	 = 18)
	plt.show()


