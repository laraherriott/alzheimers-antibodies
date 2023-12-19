# import necessary packages
import random
import matplotlib.pyplot as plt

# import necessary local classes and functions
from QSP_models.no_ab_model import NoAbModel
from QSP_models.solution import Solution

# set seed fpr reproducibility]
random.seed(1)

model = NoAbModel()

# solve the model from 0 to 1 y in time steps of 1 h
solver = Solution(model, 0, int(24*364), 1) 
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


fig1 = plt.figure(1)
plt.plot(time, brain_monomer, label = 'mon')
plt.plot(time, brain_oligomer, label = 'olig')
plt.plot(time, brain_plaque, label = 'plaque')
plt.xlabel("Time")
plt.ylabel("Species, nM")
plt.legend()
fig1.savefig('output/no_ab/change_in_brain_pre_fit.png')


fig2 = plt.figure(2)
plt.plot(time, plasma_monomer, label = 'mon')
plt.plot(time, plasma_oligomer, label = 'olig')
plt.xlabel("Time")
plt.ylabel("Species, nM")
plt.legend()
fig2.savefig('output/no_ab/change_in_plasma_pre_fit.png')

fig3 = plt.figure(3)
plt.plot(time, csf_monomer, label = 'mon')
plt.plot(time, csf_oligomer, label = 'olig')
plt.xlabel("Time")
plt.ylabel("Species, nM")
plt.legend()
fig3.savefig('output/no_ab/change_in_csf_pre_fit.png')