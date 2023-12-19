# import necessary packages
import random
import matplotlib.pyplot as plt

# import necessary local classes and functions
from QSP_models.no_ab_model import NoAbModel
from QSP_models.solution import Solution

# set seed fpr reproducibility]
random.seed(1)

model = NoAbModel()

# solve the model from 0 to 18h in time steps of 1 h
solver = Solution(model, 0, int(18), 1) 
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

total_csf = [csf_monomer[i]+csf_oligomer[i] for i in range(len(csf_oligomer))]

time_SILK = [2.116364728, 3.213552197, 4.235431889, 5.287398659, 6.377019317, 7.405204684, 8.493564208, 9.557782004, 
             10.57029326, 11.65991392, 13.77519768, 15.88579722, 16.96568912, 17.99891903, 17.98191736]
SILK_prop = [0.030717018, 0.045998002, 0.062467698, 0.077731454, 0.136426283, 0.170382111, 0.243065844, 0.301261069, 
             0.385211729, 0.454157048, 0.67332805, 0.817248389, 0.814853737, 0.876253316, 0.064880415]
SILK_values = [x*(53.2e-3+2.22e-3) for x in SILK_prop]

fig1 = plt.figure(1)
plt.plot(time, total_csf, marker='o', label='predicted')
plt.plot(time_SILK, SILK_values, marker='x', label='observed')
plt.xlabel("Time")
plt.ylabel("CSF Amyloid")
plt.legend()
fig1.savefig('output/no_ab/CSF_SILK_pre_fit.png')

time_clear = [18.9357765, 
             19.93429975, 20.96302101, 21.97167053, 22.95970586, 23.96510051, 24.96922938, 25.90301679, 26.96012778, 
             28.95988668, 30.00397818, 30.99165185, 32.00789609, 32.98743257, 33.96262922, 34.97887345, 36]

times_use = [x - 18.9357765 for x in time_clear]

prop_clear = [0.064009217, 
             0.066034605, 0.064329159, 0.064027383, 0.059395479, 0.053507292, 0.051872445, 0.048757006, 0.048101334, 
             0.040080841, 0.039337209, 0.037024124, 0.037704156, 0.033782345, 0.032231266, 0.029551349, 0.028935551]

# we now test the clearance dynamics and so turn off production and set the species concentrations to the desired steady state levels
model.params.k_in = 0
# solve the model from 0 to 10 h in time steps of 1 h
solver2 = Solution(model, 0, int(10), 1) 
solver2.y0 = [1.13, 64.72, 1300, 5.7e-3, 0.952e-3, 53.2e-3, 2.22e-3]
# will solve with LSODA method
solutions2 = solver2.solve()

time2 = solutions2.t

csf_monomer2 = solutions2.y[5]
csf_oligomer2 = solutions2.y[6]

total_csf2 = [csf_monomer2[i]+csf_oligomer2[i] for i in range(len(csf_oligomer2))]

fig2 = plt.figure(2)
plt.plot(time2, total_csf2, label='predicted')
plt.hlines(y=total_csf2[0]/2, xmin=0, xmax=10, color='k', linestyle='dashed', label='half initial concentration')
plt.xlabel("Time")
plt.ylabel("CSF Amyloid")
plt.legend()
fig2.savefig('output/no_ab/CSF_SILK_decay_pre_fit.png')
