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

# parameter values
#
# These are grouped here to simplify changing flame conditions
p          =   OneAtm               # pressure
tin_f      =   294.0                # fuel inlet temperature
tin_o      =   300.0                # oxidizer inlet temperature
#phi	   =   0.5		    # Equivalence Ratio
mdot_o     =   0.084           	    # kg/m^2/s
rho_o	   =   1.177		    # Kg/m^3 @ 300K
mdot_f     =   0.084                # kg/m^2/s
rho_f      =   0.657                # Kg/m^3 @ 294K
comp_o     =  'O2:0.21, N2:0.78, AR:0.01';   # air composition
comp_f     =  'CH4:1';                      # fuel composition

flange_length = 0.02 

alpha_0=((mdot_o/rho_o)+(mdot_f/rho_f))/(flange_length)
alpha_o=((mdot_o/rho_o)*(1+((mdot_f/mdot_o)*((rho_o/rho_f)**0.5))))/flange_length
print('Alpha_0:'+str(alpha_0))
print('Alpha_ox:'+str(alpha_o)) 
# distance between inlets is 2 cm; start with an evenly-spaced 50-point
# grid

grid_iterations=[35]#,50,100,200]

tol_ss    = [1.0e-5, 1.0e-6]        # [rtol, atol] for steady-state
                                    # problem
tol_ts    = [1.0e-5, 1.0e-2]        # [rtol, atol] for time stepping

loglevel  = 0                       # amount of diagnostic output (0
                                    # to 5)				    
refine_grid = 1                     # 1 to enable refinement, 0 to disable.

fig, (ax1,ax2,ax3)=plt.subplots(3)
ax1.set_title("Temperature Profile")
ax2.set_title("Velocity Profile")
ax3.set_title("Concentration Profile")

#h=GRI30()
#h.set(T=tin_f,P=p,X=comp_o+comp_f) #'CH4:0.5, O2:0.105,N2:0.39,AR:0.005')
#h.equilibrate('HP')
#print(h.temperature())

comptime=[]

for i in range(len(grid_iterations)):

	#initial_grid = np.array([0,0.0025,0.005,0.0075,0.008,0.009,0.01,0.0105,0.011,0.0115,0.0120,0.0125,0.01275,0.0130,0.01320,0.01340,0.01360,0.01380,0.0140,0.0141,0.0142,0.0143,0.014,0.0145,0.0146,0.0147,0.0148,0.0149,0.015,0.01525,0.01575,0.016,0.0165,0.017,0.018,0.019,0.020])
	initial_grid=np.linspace(0,1,num=grid_iterations[i])*flange_length
	dx=flange_length/len(initial_grid)
                                   # disable 				   
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
	f.solve(1)
	#elapsed=time.time()-start
	#comptime.append([grid_iterations[i],elapsed])
# Save the solution
	f.save('../results/npflame1.xml')
	
# write the velocity, temperature, and mole fractions to a CSV file
	z = f.flame.grid()
	T = f.T()
	u = f.u()
	V = f.V()
	
	umean=(np.sum(u)/len(u))	
	strain_rate=umean/dx
	print('mean strain rate'+str(strain_rate))

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

	f.showSolution()
	f.showStats()
	ax2.plot(z,u) #label='n='+str(grid_iterations[i]))
	ax1.plot(z,T)
	ax3.plot(z,V)
	#ax3.plot(z,gas.moleFraction(13), label='CH4')
	#ax3.plot(z,gas.moleFractions(4), label='O2')
	#u-profile.plot(u,T,label='n='+str(grid_iterations[i]))

#print(f.componentNames())
##########################
ax1.set_xlim(0.000, 0.020)
ax1.set_ylim(0,2500)
#ax1.set_xlabel('Z(m)')
ax1.set_ylabel('T(K)')
ax2.set_xlim(0.000, 0.020)
ax2.set_ylabel('u(m/s)')
ax3.set_ylabel('V(1/s)')
ax3.set_xlim(0.000, 0.020)
ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
plt.savefig('../plots/zvU.png')
###########################

