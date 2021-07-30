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
	
	plot_3D_angles(x, y)	
	plot_scater_heatMap('./First_order_tangent_plot.dat', 2, x, y)
	plot_scater_heatMap('./First_order_tangent_plot.dat', 3, x, y)


def plot_3D_angles(x, y):

	z = []
	with open('./Calculated_values.dat') as csv_file:
	    csv_reader = csv.reader(csv_file, delimiter=' ')
	    for row in csv_reader:
	    	z.append(float(row[0]))

	z = np.array(z)
	
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_trisurf(x, y, z, antialiased=True)
	for i in range(0,91,45):
		ax.view_init(30, i)
		plt.savefig("plot_surface_"+str(i)+".png", dpi = 100, format = "png")




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

plot_2D_func ()

