import Cantera as ct
gas=ct.GRI30('Mix')
gas.TPX=300,ct.OneAtm,'CH4:1'
fuel='CH4'

print(gas.moleFraction(str(fuel)))
#print(gas.speciesIndex('CH4'))
