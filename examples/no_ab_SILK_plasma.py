# import necessary packages
import random
import matplotlib.pyplot as plt

# import necessary local classes and functions
from QSP_models.no_ab_model import NoAbModel
from QSP_models.solution import Solution

# set seed fpr reproducibility]
random.seed(1)

model = NoAbModel()

# solve the model from 0 to 3 h in time steps of 1 h
solver = Solution(model, 0, int(3), 1)
# since we are initially testing production rates we start from 0 concentrations of each species
solver.y0 = [0, 0, 0, 0, 0, 0, 0]
# will solve with LSODA method
solutions = solver.solve()

time = solutions.t

brain_monomer = solutions.y[0]
brain_oligomer = solutions.y[1]
brain_plaque = solutions.y[2]

print(brain_monomer[-1], brain_oligomer[-1], brain_plaque[-1])

plasma_monomer = solutions.y[3]
plasma_oligomer = solutions.y[4]

print(plasma_monomer[-1], plasma_oligomer[-1])

csf_monomer = solutions.y[5]
csf_oligomer = solutions.y[6]

print(csf_monomer[-1], csf_oligomer[-1])

total_plasma = [plasma_monomer[i]+plasma_oligomer[i] for i in range(len(plasma_oligomer))]

time_SILK = [0.495211649, 1.052751315, 1.545654519, 2.04829334, 2.541999481, 3.494885456]
SILK_prop = [0.00048796, 0.025073873, 0.063703058, 0.077995403, 0.092837977, 0.091926591]
SILK_vals = [x*(5.7e-3+0.952e-3) for x in SILK_prop]

fig1 = plt.figure(1)
plt.plot(time, total_plasma, label='predicted')
plt.plot(time_SILK, SILK_vals, label='observed')
plt.xlabel("Time")
plt.ylabel("Plasma Amyloid")
plt.legend()
fig1.savefig('output/no_ab/plasma_SILK_pre_fit.png')


clear_times = [4.042087303, 5.02438086, 7.035638712, 9.032544057, 11.02503325, 13.02264116, 17.00842248, 20.99018911, 23.97862179]
clear_props = [0.081406591, 0.062295787, 0.044789249, 0.024756869, 0.020850929, 0.014290081, 0.009926749, 0.006876719, 0.005632479]

times_use = [x - 4.042087303 for x in clear_times]
vals_use = [x*(5.7e-3+0.952e-3) for x in clear_props]

# we now test the clearance dynamics and so turn off production and set the species concentrations to the desired steady state levels
model.params.k_in = 0
# solve the model from 0 to 24 h in time steps of 1 h
solver2 = Solution(model, 0, int(24), 1) 
solver2.y0 = [1.13, 64.72, 1300, 5.7e-3, 0.952e-3, 53.2e-3, 2.22e-3]
# will solve with LSODA method
solutions2 = solver2.solve()

time2 = solutions2.t

plasma_monomer2 = solutions2.y[3]
plasma_oligomer2 = solutions2.y[4]

total_plasma2 = [plasma_monomer2[i]+plasma_oligomer2[i] for i in range(len(plasma_oligomer2))]

fig2 = plt.figure(2)
plt.plot(time2, total_plasma2, label='prediction')
plt.hlines(y=total_plasma2[0]/2, xmin=0, xmax=24, color='k', linestyles='dashed', label='half initial concentration')
plt.xlabel("Time")
plt.ylabel("Plasma Amyloid")
plt.legend()
fig2.savefig('output/no_ab/plasma_SILK_decay_pre_fit.png')