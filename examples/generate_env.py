# Libraries imported to carry out the functions
import cmath
import numpy as np
from inspect import signature
from scipy.interpolate import interp1d
np.seterr(divide='ignore', invalid='ignore')

def environmentPair(L, envSeed):

	sig = signature(environmentPair)
	params = sig.parameters
	nargin = len(params)

	# These conditional statements are the same used in MATLAB
	if nargin == 1:
		env = L
	else:
		np.random.seed(envSeed)
		env = normalize(np.random.exponential(scale=1, size=(2, L)))

	lenv = np.zeros(shape = env.shape, dtype=complex)
	R,C = lenv.shape
	for r in range(R):
		for c in range(C):
			lenv[r,c] = cmath.log(env[r,c])

	a = lenv[0, :]
	b = lenv[1, :]

	x = np.arange(0, 5, 0.01)

	x2env = lambda x:normalize(np.exp(((a+b)/2+x/2*np.array([[a-b], [b-a]]))))
	env2dE = lambda e:np.linalg.norm(e[0,:]-e[1,:])
	x2dE = [env2dE(x2env(val)) for val in x]

	lim = np.where(np.diff(x2dE)==0)
	if not len(lim[0]) == 0:
		x2dE = np.delete(x2dE, range(lim[0][0], len(x2dE)))
		x = np.delete(x, range(lim[0][0], len(x)))

	envPair = lambda dE:x2env(interp1d(x2dE, x)(dE))

	return envPair

def normalize(arr):
	if len(arr.shape) > 2:
		arr = arr[:, 0, :]
	temp = arr@arr.T
	diags = np.diag(temp)
	diags = diags[:, np.newaxis]
	normalized = arr/np.sqrt(diags)
	return normalized


