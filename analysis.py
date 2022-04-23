from result_analysis.plot_UEF import *

print("ok")
r_exp, r_var = get_UEF_data(1) 
print("plot")
plot_UEfrontier(r_exp, r_var)