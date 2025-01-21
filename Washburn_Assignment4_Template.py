# Script modified by Tucker J Lancaster
import argparse, random
import numpy as np
import matplotlib.pyplot as plt

"""
Assignment 4: Bacterial Warfare 

General instructions: Please use this code skeleton for your implementation, making no changes to the function names,
function parameter names, or function return-value names. We recommend that you complete this assignment in the
following order:
	1. implement the argument parser, as described in the "if __name__ == '__main__'" block at the bottom of this script
	2. implement seed_grid
	3. implement update_grid
	4. implement update_plot
	5. run your code and confirm that the plots you see resemble those from the paper
"""

def seed_grid(n): 
	"""
	This function is used to create the initial "playing field", an n x n array randomly populated with 1's (for the
	red species) and -1's (for the blue species). Initially, populate every position in the grid (i.e., fill every
	position with either a 1 or a -1).

	:param n: grid dimension. returned grid will be nxn
	:type n: int
	:return grid: n x n grid, randomly populated with 1's and -1's
	:rtype grid: numpy.ndarray
	"""
	
	#https://www.w3schools.com/python/numpy/numpy_random.asp#:~:text=The%20choice()%20method%20allows,returns%20one%20of%20the%20values.

	grid = np.random.choice([-1,1], (n,n)).astype(np.int8)
	
	return grid

def update_grid(grid, killing_rate, reproducing_rate):
	"""
	This function is used to simulate a single step of the bacterial warfare simulation. It takes in the current
	playing field (grid), updates it (through killing and reproduction), and returns the updated playing field (grid).

	:param grid: initial grid
	:type grid: numpy.ndarray
	:param killing_rate: proportion of cells that activate their TS66 system each step
	:type killing_rate: float
	:param reproducing_rate: proportion of cells that attempt to reproduce each step
	:type reproducing_rate: float
	:return grid: updated grid
	:rtype grid: numpy.ndarray
	"""

	

	# Determine how many bacteria are activated for T6SS (killing_rate * number of cells alive)
	# alive = np.count_nonzero(grid == 1) + np.count_nonzero(grid == -1)
	alive = np.count_nonzero(grid)
	num_activated = int(killing_rate * alive) # Properly define this variable 
	

	# use a for loop to activate the specified number of bacteria.
	for i in range(num_activated):
		# Randomly choose a grid position
		x = random.randint(0, (len(grid)-1))  # define x appropriately
		y = random.randint(0, (len(grid)-1))  # define y appropriately
		# using a while loop, choose new random grid positions until you encounter an occupied position (if the position
		# you chose initially was unoccupied)
		while grid[x, y] == 0:
			x = random.randint(0, (len(grid)-1))  # define x appropriately 
			y = random.randint(0, (len(grid)-1))  # define y appropriately 

		
		
		# kill all bacteria of the opposite species occupying positions adjacent to the bacterium in the chosen position
		#add code here to do this - will set the opposite species to zero
		#maybe slice grid to only select the area around the (x,y)
		#if cell is not the same as (x,y) then kill it
		#https://www.digitalocean.com/community/tutorials/python-numpy-where
			#might be useful to determine if values around it are the same or not and what to replace it as 
		if (x-1) > 0:
			low_x = x-1
		else:
			low_x = 0

		if (x+1) > len(grid)-1:
			high_x = len(grid)
		else:
			high_x = x+2

		if (y-1) > 0:
			low_y = y-1
		else:
			low_y = 0
		
		if (y+1) > len(grid)-1:
			high_y = len(grid)
		else:
			high_y = y+2

		adjacents = grid[low_x:high_x, low_y:high_y]
		adjacents[adjacents == (grid[x,y] * -1)] = 0
		 



	# Determine the number of bacteria that will reproduce (number of bacteria alive * reproduction rate)
	#https://stackoverflow.com/questions/28663856/how-do-i-count-the-occurrence-of-a-certain-item-in-an-ndarray
	alive = np.count_nonzero(grid == 1) + np.count_nonzero(grid == -1) #number of alive bacteria
	num_reproducing = int(alive * reproducing_rate) # Properly define this variable 

	# Loop through number of bacteria that are reproducing
	for i in range(num_reproducing):
		# Randomly choose a grid position
		x = random.randint(0, (len(grid)-1))  # define x appropriately
		y = random.randint(0, (len(grid)-1))  # define y appropriately
		# using a while loop, choose new random grid positions until you encounter an occupied position (if the position
		# you chose initially was unoccupied)
		while grid[x, y] == 0:
			x = random.randint(0, (len(grid)-1))  # define x appropriately
			y = random.randint(0, (len(grid)-1))  # define y appropriately

		

		# if there is not an empty position adjacent to the cell, continue without modifying the grid
		

		if (x-1) > 0:
			low_x = x-1
		else:
			low_x = 0

		if (x+1) > len(grid)-1:
			high_x = len(grid)
		else:
			high_x = x+2

		if (y-1) > 0:
			low_y = y-1
		else:
			low_y = 0
		
		if (y+1) > len(grid)-1:
			high_y = len(grid)
		else:
			high_y = y+2

		reproduce = grid[low_x:high_x, low_y:high_y]
		
		#https://stackoverflow.com/questions/73898808/how-to-randomly-edit-elements-in-a-2d-array
		if np.any(reproduce == 0):
			#randomly reproduce
			if grid[x,y] == 1 or grid[x,y] == -1:
				#replace 0 with 1 or -1
				j,k = np.where(reproduce == 0)
				random_select = np.random.choice(np.arange(len(j)), 1)
				reproduce[j[random_select], k[random_select]] = grid[x,y]
				
			
			
		# if there is at least one empty position adjacent to the cell, set one of the empty positions to the species
		# of the reproducing cell. The filled position should be chosen randomly from all available positions adjacent
		# to the focal bacterium.

	

	return grid

def update_plot(img_obj, grid, i):
	"""
	This function is used to update the plot of the current bacterial distribution across the playing field. The
	parameter img_obj should be a matplotlib.image.AxesImage object, which is the type of object returned by
	matplotlib.pyplot.matshow().

	:param img_obj: plotting object to be updated
	:type img_obj: matplotlib.image.AxesImage
	:param grid: array containing the current playing field
	:type grid: numpy.ndarray
	:param i: current simulation step
	:type i: int
	"""

	# use the matplotlib.image.AxesImage.set_data() method to update array associated with img_obj
	plt.matplotlib.image.AxesImage.set_data(img_obj, grid)
	# update the title to be "Step i", where i is the current step
	plt.title(f"Step {i}")

	# use the matplotlib.pyplot.draw() command to redraw the figure
	plt.matplotlib.pyplot.draw()
	# use the matplotlib.pyplot.pause() command to delay the drawing of the next frame by 0.1 seconds
	plt.matplotlib.pyplot.pause(0.1)
	  # remove this pass statement when your function is implemented

	# this return statement is not necessary, but is included to help with grading
	return img_obj


if __name__ == '__main__':
	# Use argparse to take in the grid size, the number steps you want to simulate, the killing_rate, and the
	# reproduction_rate. Set the default values of these arguments to what was used in the paper. Use flagged arguments
	# for all four arguments. Use the flag names "--gridsize", "--n_steps", "--killing_rate", and "--reproduction_rate".
	# include a descriptive help string for each, and specify the expected type. Please place all of your argument
	# parsing code here, within the __name__ == '__main__' block.

	# add your argument parsing code here
	parser = argparse.ArgumentParser(description = 'Using python to study the effect of killing on spatial structure in bacterial populations.') #fix this
	parser.add_argument('--gridsize', help='the size of grid', default = 500, type=int)
	parser.add_argument('--n_steps', help = 'the number of steps to simulate', default = 100, type = int) 
	parser.add_argument('--killing_rate', help = 'killing rate of bacteria', default = 0.05, type = float)
	parser.add_argument('--reproduction_rate', help = 'reproduction rate of bacteria', default = 0.05, type = float)
	args = parser.parse_args()  # properly define this variable

	# the code below is complete, but you should read it to better understand how the functions you write above should
	# operate
	grid = seed_grid(args.gridsize)
	img_obj = plt.matshow(grid, cmap='seismic')
	for i in range(args.n_steps):
		print('Loop: ' + str(i))
		update_grid(grid, args.killing_rate, args.reproduction_rate)
		update_plot(img_obj, grid, i)





