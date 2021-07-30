import csv
import numpy as np
import matplotlib.pyplot as plt

def plot_2D_func ():

	x = []
	y = []


	with open('./Mesh_points.dat') as csv_file:
		csv_reader = csv.reader(csv_file, delimiter=' ')
		for row in csv_reader:
		    x.append(float(row[0]))
		    y.append(float(row[1]))

	x = np.array(x)
	y = np.array(y)
	
	report = open("./report.tex", "w")
	tex_main = open("./src/tex/tex_body","r") 
	tex_body = tex_main.readlines()
	tex_main.close()
	
	tex_inc_plot = open("./src/tex/tex_inc_plot","r") 
	tex_plot = tex_inc_plot.readlines()
	tex_inc_plot.close()
	
	report.writelines(tex_body[:6])	

	name_of_plots = plot_3D_angles(x, y)
	
	for i in name_of_plots:
		report.write("\n") 
		report.writelines(tex_plot[:2])
		report.write("\includegraphics[width=0.9\\textwidth]{./" + i + "}\n") 
		report.write("\caption{surface plot.}\n")
		report.writelines(tex_plot[4:])
		
	
	
	for i in range(2,4):	
	
		name_of_plot = plot_scater_heatMap('./First_order_tangent_plot.dat', i, x, y)	

		report.write("\n") 
		report.writelines(tex_plot[:2])
		report.write("\includegraphics[width=0.9\\textwidth]{./" + name_of_plot + "}\n") 
		report.write("\caption{Heat map for $x_" + str(i-2) + "$.}\n")
		report.writelines(tex_plot[4:])
		
	
	report.writelines(tex_body[6:])
	report.close()


def plot_3D_angles(x, y):

	z = []
	with open('./Calculated_values.dat') as csv_file:
	    csv_reader = csv.reader(csv_file, delimiter=' ')
	    for row in csv_reader:
	    	z.append(float(row[0]))

	z = np.array(z)
	
	
	name_of_plots = []
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_trisurf(x, y, z, antialiased=True)
	for i in range(0,91,45):
		ax.view_init(30, i)
		plt.savefig("plot_surface_"+str(i)+".png", dpi = 100, format = "png")
		name_of_plots.append("plot_surface_"+str(i)+".png")
	
	plt.show()
	return name_of_plots



def plot_scater_heatMap(fileName, col, x, y):
	DPI = 100
	z = []
	with open(fileName) as csv_file:
	    csv_reader = csv.reader(csv_file, delimiter=' ')
	    for row in csv_reader:
	    	z.append(float(row[col]))

	z = np.array(z)
		
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.scatter(x, y, s=(2*DPI**2/x.size), c=z)
	plt.savefig("plot_scater_heatMap_"+str(col)+".png", dpi = DPI, format = "png")

	return "plot_scater_heatMap_"+str(col)+".png"

plot_2D_func ()

