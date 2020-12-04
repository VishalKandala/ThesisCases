#!usr/bin/python
from Cantera import *
from Cantera.OneD import *
from Cantera.OneD.CounterFlame import CounterFlame
from Cantera.num import array
import numpy as np
import matplotlib.pyplot as plt
import time
tol_ss    = [1.0e-5, 1.0e-6]        # [rtol, atol] f1or steady-state
                                    # problem
tol_ts    = [1.0e-5, 1.0e-2] 

d=0.02

# Flow Parameters #
p          =   OneAtm               # pressure
tin_f      =   294.0                # fuel inlet temperature
tin_o      =   300.0                # oxidizer inlet temperature
phi	   =   3.17         # Equivalence Ratio     
rho_o	   =   1.177		    # Kg/m^3 @ 300K
mdot_f     =   0.084                # kg/m^2/s
rho_f      =   (0.657*(phi/(2+phi)))+(rho_o*2/(2+phi)) # Kg/m^3 @ 294K
mdot_o     =  0.084                 # Kg/m^2/s
aircomp=[0.21,0.78,0.01]
comp_o     =  'O2:0.21, N2:0.78, AR:0.01';  # air composition
comp_f     =  'CH4:'+str(phi/(2+phi))+', O2:'+str(2*aircomp[0]/(2+phi))+', N2:'+str(2*aircomp[1]/(2+phi))+', AR:'+str(2*aircomp[2]/(2+phi))
initial_grid=d*np.linspace(0,1,num=10)

gas = GRI30('Mix')
gas.addTransportModel('Multi')

# create an object representing the counterf1low flame configuration,
# which consists of1 a fuel inlet on the left, the flow in the middle,
# and the oxidizer inlet on the right. Class CounterFlame creates this
# conf1iguration.
f1 = CounterFlame(gas = gas, grid = initial_grid)
f1.fuel_inlet.set(massflux = mdot_f,mole_fractions = comp_f,temperature = tin_f)
f1.oxidizer_inlet.set(massflux = mdot_o,mole_fractions = comp_o,temperature = tin_o)
# set the error tolerances
f1.set(tol = tol_ss, tol_time = tol_ts)
# construct the initial solution estimate. To do so, it is necessary
# to specif1y the fuel species. If a fuel mixture is being used,
# specif1y a representative species here for the purpose of
# constructing an initial guess.
f1.init(fuel = 'CH4')
filename='../results/CH4_Air_a_3.17_energy_multi_soret.xml'
f1.restore(file=filename, id='solution')

f2 = CounterFlame(gas = gas, grid = initial_grid)
f2.fuel_inlet.set(massflux = mdot_f,mole_fractions = comp_f,temperature = tin_f)
f2.oxidizer_inlet.set(massflux = mdot_o,mole_fractions = comp_o,temperature = tin_o)
# set the error tolerances
f2.set(tol = tol_ss, tol_time = tol_ts)
# construct the initial solution estimate. To do so, it is necessary
# to specif1y the fuel species. If a fuel mixture is being used,
# specif1y a representative species here for the purpose of
# constructing an initial guess.
f2.init(fuel = 'CH4')
filename='../results/CH4_Air_a_3.17_energy_multi.xml'
f2.restore(file=filename, id='solution')

plt.plot(f1.flame.grid(),f1.T())
plt.plot(f2.flame.grid(),f2.T())
plt.savefig('../plots/tempcomp.png')
