import csv
import numpy as np
import matplotlib.pyplot as plt

def plot_1D_func ():

	x = []


	with open('./Mesh_points.dat') as csv_file:
		csv_reader = csv.reader(csv_file, delimiter=' ')
		for row in csv_reader:
		    x.append(float(row[0]))

	x = np.array(x)
	
	
	y = []
	with open('./Calculated_values.dat') as csv_file:
	    csv_reader = csv.reader(csv_file, delimiter=' ')
	    for row in csv_reader:
	    	y.append(float(row[0]))

	y = np.array(y)
	
	report = open("./report.tex", "w")
	tex_main = open("./src/tex/tex_body","r") 
	tex_body = tex_main.readlines()
	tex_main.close()
	
	tex_inc_plot = open("./src/tex/tex_inc_plot","r") 
	tex_plot = tex_inc_plot.readlines()
	tex_inc_plot.close()
	
	report.writelines(tex_body[:6])	

	
	
	name_of_plot = plot_scater_heatMap('./First_order_tangent_plot.dat', 1, x, y)	

	report.write("\n") 
	report.writelines(tex_plot[:2])
	report.write("\includegraphics[width=0.9\\textwidth]{./" + name_of_plot + "}\n") 
	report.write("\caption{Grad map for $x_0$.}\n")
	report.writelines(tex_plot[4:])
		
	
	report.writelines(tex_body[6:])
	report.close()



def plot_scater_heatMap(fileName, col, x, y):
	_scale = 200000
	z = []
	with open(fileName) as csv_file:
	    csv_reader = csv.reader(csv_file, delimiter=' ')
	    for row in csv_reader:
	    	z.append(float(row[col]))

	z = np.array(z)
		
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.scatter(x, y, s=(_scale/x.size), c=z)
	plt.savefig("plot_scater_heatMap_"+str(col)+".png", dpi = 100, format = "png")

	return "plot_scater_heatMap_"+str(col)+".png"

plot_2D_func ()

