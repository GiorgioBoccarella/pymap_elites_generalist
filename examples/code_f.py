# Libraries imported to carry out the functions
import cmath # Used to caary out complex math
import numpy as np # Used to manipulate matrices
from inspect import signature # Used to find the number of arguments of a function
from scipy.interpolate import interp1d # Used to interpolate 1D data
from scipy.stats import zscore as normalize # Used to normalize
np.seterr(divide='ignore', invalid='ignore')
def environmentPair(L, envSeed):
	'''
	Definition of the funciton that will be used.
	'''

	# The follwing lines are used to calculate the number of arguments that the function has. The number 
	# of arguments in MATLAB can are stored in a variable called "nargin" but there is no such variable in 
	# python therefore we use 3rd party libraries to find the number of parameters which are then stored in a
	# variable named "nargin".
	sig = signature(environmentPair)
	params = sig.parameters
	nargin = len(params)

	# These conditional statements are the same the ones used in MATLAB
	if nargin == 1:
		env = L
	else:
		np.random.seed(envSeed) # In MATLAB "rng" is used to make sure the same random number are used everytime, in python we use "seed"
		# The same funcitons are applied here as used in MATLAB but in python these functions are provided by 3rd part libraries
		env = normalize(np.random.exponential(scale = 1, size = (2, L)), ddof = 1)
	# In the following lines we declare a complex array of zeros and then we store the logs of variable "env" in the variable of zeros
	# In order to calculate the log with complex numbers we use cmath library. 
	lenv = np.zeros(shape = env.shape, dtype = complex)
	R,C = lenv.shape
	for r in range(R):
		for c in range(C):
			lenv[r,c] = cmath.log(env[r,c])

	# This is similar to the MATLAB code but in python indices start from 0 however in MATLAB they start from 1.
	a = lenv[0,:]
	b = lenv[1,:]

	# In this part a range of values are stored in a variables like we saw in MATLAB.
	x = np.arange(0, 5, 0.01)

	# In MATLAB function handles are defined by "@" sign, however, in python we define function handles using "lambda". Using
	# lambda we carry out the same tasks as we saw in the MATLAB code. The mechanics are the same, the only difference is in
	# the function nomencletrue and base libraries.
	x2env = lambda x:normalize(np.exp(((a+b)/2+x/2*np.array([[a-b], [b-a]])))) # Similar to the matlab function but written in ptython syntax
	env2dE = lambda e:quartNorm(e[0,:]-e[1,:]) # There is no quarternion Norm funciton in python therefore we defined our own.
	x2dE = [env2dE(x2env(val)) for val in x] # There is no "arrayfun" function in python, therefore we use for loops to carry out this task.

	# In the following lines first we find the difference between adjacent elements of the array and then we find which of the elements in 
	# difference array element is equal to zero. If there is an element that does have a value of zero then we delete all the elements in "x"
	# and x2dE variables from the first occurance of zero till the end. 
	lim = np.where(np.diff(x2dE)==0) #The location of zeros are calculated in the difference array
	if not len(lim[0]) == 0:
		 # All the values from the first occurance of zero till the end are deleted from both of the following arrays
		x2dE = np.delete(x2dE, range(lim[0][0], len(x2dE)))
		x = np.delete(x, range(lim[0][0], len(x)))

	# As seen in the original MATLAB function, another function handle is created that is then returned to the output.
	envPair = lambda dE:x2env(interp1d(x2dE, x)(dE))

	return envPair

def quartNorm(arr):
	'''
	This funciton gets the quarternion norm for the passed array. This is based on the MATLAB "norm" function. For more detail look at the
	documentation of the MATLAB "norm" function.
	'''
	R, C = arr.shape
	SUM = 0
	for c in range(C):
		SUM += arr[0,c]**2
	return np.sqrt(SUM)


for i in range(0, 5):
	for j in range(0, 5):
		print(i, j)
		envPair = environmentPair(i, j)
		print(envPair)


envPair(1.3)