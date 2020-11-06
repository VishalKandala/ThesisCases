#!/usr/bin/python
# NPFLAME1 - A nonpremixed counterflow flame.
#
#    This script computes an atmospheric-pressure methane/air
#    counterflow flame using GRI-Mech 3.0.

from Cantera import *
from Cantera.OneD import *
from Cantera.OneD.CounterFlame import CounterFlame
from Cantera.num import array
import numpy as np
import matplotlib.pyplot as plt
import time
####### Solution Controls ######################################

tol_ss    = [1.0e-5, 1.0e-6]        # [rtol, atol] for steady-state
                                    # problem
tol_ts    = [1.0e-5, 1.0e-2]        # [rtol, atol] for time stepping

loglevel  = 0                       # amount of diagnostic output (0
                                    # to 5)				    
refine_grid = 1                     # 1 to enable refinement, 0 to disable.

###### Grid Initialization ######################################
d1=2 #flange distance in cm.
d=d1*0.01 # flange distance in m.

phi_range=[1.8,2.2,3.17]
grid_size=[10,20,30]

# Iterating over  various equivalence ratios
for i in range(len(phi_range)):
############ Flame Parameters ###############################
	p          =   OneAtm               # pressure
	tin_f      =   294.0                # fuel inlet temperature
	tin_o      =   300.0                # oxidizer inlet temperature
	phi	   =   phi_range[i]         # Equivalence Ratio     
	rho_o	   =   1.177		    # Kg/m^3 @ 300K
	mdot_f     =   0.084                # kg/m^2/s
	rho_f      =   (0.657*(phi/(2+phi)))+(rho_o*2/(2+phi)) # Kg/m^3 @ 294K
	mdot_o     =  0.084                 # Kg/m^2/s
	aircomp=[0.21,0.78,0.01]
	comp_o     =  'O2:0.21, N2:0.78, AR:0.01';  # air composition
	comp_f     =  'CH4:'+str(phi/(2+phi))+', O2:'+str(2*aircomp[0]/(2+phi))+', N2:'+str(2*aircomp[1]/(2+phi))+', AR:'+str(2*aircomp[2]/(2+phi))
	initial_grid=np.linspace(0,1,num=grid_size[i])*d 
	# Cold Strain Rates # 

	alpha_0=((mdot_o/rho_o)+(mdot_f/rho_f))/(d)
	alpha_o=((mdot_o/rho_o)*(1+((mdot_f/mdot_o)*((rho_o/rho_f)**0.5))))/d			   
	print('Phi:'+str(phi))
	print('Alpha_0:'+str(alpha_0))
	print('Alpha_ox:'+str(alpha_o))

################ create the gas object ########################
#
# This object will be used to evaluate all thermodynamic, kinetic,
# and transport properties
#

# Here we use GRI-Mech 3.0 with mixture-averaged transport
# properties. To use your own mechanism, use function
# IdealGasMix('mech.cti') to read a mechanism in Cantera format.  If
# you need to convert from Chemkin format, use the ck2cti utility
# program first.
	gas = GRI30('Mix')

# create an object representing the counterflow flame configuration,
# which consists of a fuel inlet on the left, the flow in the middle,
# and the oxidizer inlet on the right. Class CounterFlame creates this
# configuration.

	f = CounterFlame(gas = gas, grid = initial_grid)

# Set the state of the two inlets

	f.fuel_inlet.set(massflux = mdot_f,
                 mole_fractions = comp_f,
                 temperature = tin_f)

	f.oxidizer_inlet.set(massflux = mdot_o,
                     mole_fractions = comp_o,
                     temperature = tin_o)

# set the error tolerances
	f.set(tol = tol_ss, tol_time = tol_ts)

# construct the initial solution estimate. To do so, it is necessary
# to specify the fuel species. If a fuel mixture is being used,
# specify a representative species here for the purpose of
# constructing an initial guess.
	f.init(fuel = 'CH4')

# show the starting estimate
#	f.showSolution()

# First disable the energy equation and solve the problem without
# refining the grid
	start=time.time()
	f.setRefineCriteria(ratio = 3, slope = 1, curve = 1, prune = 0)	
	f.set(energy = 'off')
	start=time.time()	
	f.solve(loglevel, 1)
# Now specify grid refinement criteria, turn on the energy equation,
# and solve the problem again. The ratio parameter controls the
# maximum size ratio between adjacent cells; slope and curve should be
# between 0 and 1 and control adding points in regions of high
# gradients and high curvature, respectively. If prune > 0, points
# will be removed if the relative slope and curvature for all
# components fall below the prune level. Set prune < min(slope,
# curve), or to zero to disable removing grid points.
	f.setRefineCriteria(ratio = 3, slope = 0.1, curve = 0.1, prune = 0)
	f.set(energy = 'on')
	f.solve(loglevel,1)
	elapsed=time.time()-start
	print('Solution Time:'+str(elapsed))
# Save the solution
	f.save('../results/CH4_Air_'+str(phi)+'_'+str(d1).zfill(3)+'.xml')
	
# write the velocity, temperature, and mole fractions to a CSV file
	z = f.flame.grid()
	T = f.T()
	u = f.u()
	V = f.V()
	vmax=np.amax(V)

	#results=np.array([z,T,u,V]) # saving all the arrays into a single 2-D array.
	#filename='../results/npflame1_ref0_'+str(grid_iterations[i])+'.csv'
	#np.savetxt(filename, results, delimiter=',')
#	fcsv = open('../results/npflame1_'+str(grid_iterations[i])+'.csv',"w")
#	writeCSV(fcsv, ['z (m)', 'u (m/s)', 'V (1/s)', 'T (K)']+list(gas.speciesNames()))
#	for n in range(f.flame.nPoints()):
#    		f.setGasState(n)
#    		writeCSV(fcsv, [z[n], u[n], V[n], T[n]]+list(gas.moleFractions()))
#		fcsv.close()
#	print 'solution saved to :'+filename

#	f.showSolution()
#	f.showStats()
############# Plots ##########
	fig,(ax1,ax2,ax3)=plt.subplots(3)
	fig.suptitle('$\phi = '+str(phi)+', d = '+str(d1)+'cm, t = '+str(elapsed)[:4]+'s'+', '+'\alpha = '+str(alpha_o)[:3]+'(1/s), V = '+str(vmax)+'$')
	ax2.plot(z,u,color='b')
	ax1.plot(z,T,color='r')
	ax3.plot(z,V,color='g')
	ax1.set_xlim(0.000, 0.020)
	ax1.set_ylim(0,2500)
	ax1.set_ylabel('T(K)')
	ax2.set_xlim(0.000, 0.020)
	ax2.set_ylabel('u(m/s)')
	ax3.set_ylabel('V(1/s)')
	ax3.set_xlim(0.000, 0.020)
	ax1.grid(True)
	ax2.grid(True)
	ax3.grid(True)
	plt.savefig('../plots/CH4-Air_1D_'+str(phi)+'_'+str(d1).zfill(3)+'.png')
