#
# Class to solve the differential equations
#
import scipy

class Solution:
    def __init__(self, model, t_0, t_end, step_size):
        self.t_start = t_0
        self.t_end = t_end
        self.t_eval = list(range(t_0, t_end, step_size))
        self.model = model
        if self.model.type == 'no_ab':
            #self.y0 = [0, 0, 0, 0, 0, 0, 0]
            self.equations = model.brain
            self.y0 = [1.13, 64.72, 1300, 5.7e-3, 0.952e-3, 53.2e-3, 2.22e-3]
        elif self.model.type == 'one_ab':
            self.y0 = [1.13, 1.86, 1300, 5.7e-3, 0.952e-3, 53.2e-3, 2.22e-3]

    def solve(self):
        solution = scipy.integrate.solve_ivp(fun=lambda t, y: self.equations(t, y),
                                             t_span=[self.t_eval[0],
                                                     self.t_eval[-1]],
                                             y0=self.y0,
                                             t_eval=self.t_eval,
                                             max_step=0.1,
                                             method='LSODA')
        return solution
