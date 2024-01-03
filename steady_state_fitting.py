import random
import numpy as np
import scipy
import pypesto
import pypesto.optimize as optimize
import matplotlib.pyplot as plt

# import necessary local classes and functions
from QSP_models.no_ab_model import NoAbModel
from QSP_models.solution import Solution

csf_observed = [0.077731454, 0.136426283, 0.170382111, 0.243065844, 0.301261069, 0.385211729, 0.454157048] # t = 5, 6, 7, 8, 9, 10, 11
csf_silk = [x*(53.2e-3+2.22e-3) for x in csf_observed]
plasma_observed = [0.025073873, 0.077995403, 0.091926591]
plasma_silk = [x*(5.7e-3+0.952e-3) for x in plasma_observed]

csf_clear = [0.040080841, 0.039337209, 0.037024124, 0.037704156, 0.033782345, 0.032231266, 0.029551349, 0.028935551]
csf_clear_silk = [x*(53.2e-3+2.22e-3) for x in csf_clear]
plasma_clear = [0.009926749, 0.006876719, 0.005632479]
plasma_clear_silk = [x*(5.7e-3+0.952e-3) for x in plasma_clear]

def diff(a, b):
    l = list()
    for x, y in zip(a, b):
        l.append(x - y)

    return l

def sum_of_squares(l):
    squares = [x**2 for x in l]
    return sum(squares)

def relative_sum_of_squares(a, b):
    squares = []
    for x, y in zip(a, b):
        if y != 0:
            squares.append(((y-x)/x)**2)
        else:
            squares.append(1)
    return sum(squares)

def error(x):
    model_1 = NoAbModel()
    model_1.params.k_olig_inc = x[0]
    model_1.params.k_olig_sep = x[1]
    model_1.params.k_clear_Abeta_brain = x[2]
    model_1.params.k_clear_Abeta_plasma = x[3]
    model_1.params.k_clear_Abeta_csf = x[4]
    model_1.params.k_clear_oligomer_brain = x[5]
    model_1.params.k_clear_oligomer_plasma = x[6]
    model_1.params.k_clear_oligomer_csf = x[7]
    model_1.params.k_monomer_brain_plasma = x[8]
    model_1.params.k_monomer_plasma_brain = x[9]
    model_1.params.k_monomer_brain_csf = x[10]
    model_1.params.k_plaque_inc = x[11]
    model_1.params.k_plaque_sep = x[12]
    model_1.params.k_oligomer_brain_plasma = x[13]
    model_1.params.k_oligomer_plasma_brain = x[14]
    model_1.params.k_oligomer_brain_csf = x[15]
    model_1.params.k_monomer_plasma_csf = x[16]
    model_1.params.k_monomer_csf_plasma = x[17]
    model_1.params.k_oligomer_plasma_csf = x[18]
    model_1.params.k_oligomer_csf_plasma = x[19]

    # to test production and compare to SILK
    solver_1 = Solution(model_1, 0, int(24), 1)
    solver_1.y0 = [0, 0, 0, 0, 0, 0, 0]
    solutions_1 = solver_1.solve()

    csf = [solutions_1.y[5][i] + solutions_1.y[6][i] for i in range(len(solutions_1.y[5]))] # total csf conc
    csf_use = csf[4:12]
    plasma = [solutions_1.y[3][i] + solutions_1.y[4][i] for i in range(len(solutions_1.y[3]))] # total plasma conc
    plasma_use = plasma[1:4]

    # csf_diff = diff(csf_use, csf_silk)
    # plasma_diff = diff(plasma_use, plasma_silk)

    csf_relative = relative_sum_of_squares(csf_silk, csf_use)
    plasma_relative = relative_sum_of_squares(plasma_silk, plasma_use)

    # to test maintenance of steady state over a year (with slight increase in plaque)
    solver_2 = Solution(model_1, 0, int(24*364), 1)  
    solver_2.y0 = [1.13, 64.72, 1300, 5.7e-3, 0.952e-3, 53.2e-3, 2.22e-3]
    solutions_2 = solver_2.solve()

    csf_ss1 = solutions_2.y[5][-1]
    csf_ss2 = solutions_2.y[6][-1] # total csf conc
    plasma_ss1 = solutions_2.y[3][-1]
    plasma_ss2 = solutions_1.y[4][-1] # total plasma conc
    monomer_ss = solutions_2.y[0]
    oligomer_ss = solutions_2.y[1]
    plaque_ss = solutions_2.y[2]

    final_vals = [plaque_ss[-1], oligomer_ss[-1], monomer_ss[-1], plasma_ss1, plasma_ss2, csf_ss1, csf_ss2]
    desired = [1300, 64.72, 1.12, 5.7e-3, 0.952e-3, 53.2e-3, 2.22e-3]

    steady_state_error = 0
    
    for value, observed in zip(final_vals, desired):
        if value != 0:
            steady_state_error += ((value-observed/observed))**2
        else:
            steady_state_error += 1

    # to test clearance and compare to SILK
    model_1.params.k_in = 0

    solver_3 = Solution(model_1, 0, int(36), 1)
    solver_3.y0 = [1.13, 64.72, 1300, 5.7e-3, 0.952e-3, 53.2e-3, 2.22e-3]
    solutions_3 = solver_3.solve()

    csf_clear = [solutions_3.y[5][i] + solutions_3.y[6][i] for i in range(len(solutions_3.y[5]))] # total csf conc
    csf_clear_use = csf_clear[29:]
    plasma_clear = [solutions_3.y[3][i] + solutions_3.y[4][i] for i in range(len(solutions_3.y[3]))] # total plasma conc
    plasma_clear_use = [plasma_clear[17], plasma_clear[21], plasma_clear[24]]

    # csf_clear_diff = diff(csf_clear_use, csf_clear_silk)
    # plasma_clear_diff = diff(plasma_clear_use, plasma_clear_silk)

    csf_clear_relative = relative_sum_of_squares(csf_clear_silk, csf_clear_use)
    plasma_clear_relative = relative_sum_of_squares(plasma_clear_silk, plasma_clear_use)

    # overall = csf_diff + plasma_diff + csf_clear_diff + plasma_clear_diff

    return csf_relative + plasma_relative + plasma_clear_relative + csf_clear_relative + 4*steady_state_error


objective = pypesto.Objective(fun=error)#, res=error)

# parameter order:
# k_olig_inc, k_olig_sep, k_clear_Abeta_brain, k_clear_Abeta_plasma, k_clear_Abeta_csf
# k_clear_oligomer_brain, k_clear_oligomer_plasma,k_clear_oligomer_csf, k_monomer_brain_plasma
# k_monomer_plasma_brain, k_monomer_brain_csf, k_plaque_inc, k_plaque_sep
# k_oligomer_brain_plasma, k_oligomer_plasma_brain, k_oligomer_brain_csf,
# k_monomer_plasma_csf, k_monomer_csf_plasma, k_oligomer_plasma_csf, k_oligomer_csf_plasma

lb = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
ub = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

# initial estimates of parameter values
start = np.array([0.924, 0, 0.165, 0.231, 0, 0.165, 0.231, 0, 0.0923, 0, 0.0982, 0.165, 0, 0.0923, 0, 0.0982, 0, 0.065, 0, 0.065])
# start = np.array([4.11026193e-02, 9.94749664e-01, 9.64366261e-01, 1.00000000e+00,
#        2.27030315e-01, 1.00000000e+00, 9.99979441e-01, 8.78256393e-01,
#        6.54693652e-03, 9.47700481e-01, 3.30419301e-02, 0.00000000e+00,
#        1.00000000e+00, 6.68065993e-01, 1.07472400e-01, 1.69385523e-04,
#        9.93490129e-01, 2.15326288e-03, 1.54677175e-01, 9.99455278e-01])
start = start.reshape((1, 20))
problem = pypesto.Problem(objective=objective, lb=lb, ub=ub, x_guesses=start)#, startpoint_method=custom_startpoints)

optimizer = optimize.ScipyOptimizer(options = {'maxfun': 2500})

result = optimize.minimize(problem=problem, optimizer=optimizer, n_starts=10)

print(result.optimize_result.list[0])

file_object = open('no_antibody.txt', 'w')
file_object.write(str(result.optimize_result.list[0]))
file_object.close()
