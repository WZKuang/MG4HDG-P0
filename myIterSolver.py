from ngsolve import *
from ngsolve.krylovspace import LinearSolver

class IterSolver(LinearSolver):
    def __init__(self, *args, **kargs):
        super().__init__(*args, **kargs)

    def _SolveImpl(self, rhs: BaseVector, sol: BaseVector):
        r = rhs.CreateVector()
        d = rhs.CreateVector()
        r.data = rhs - self.mat * sol
        d.data = self.pre * r
        # A_tilder norm for residual checking here - A assumed SPD
        res_norm = sqrt(InnerProduct(d, r))
        if self.CheckResidual(res_norm):
            return
        
        while True:
            # self.pre => approximation of A inverse
            sol.data += d
            r.data = rhs - self.mat * sol
            d.data = self.pre * r
            prev_res_norm = res_norm
            res_norm = sqrt(InnerProduct(d, r))
            if res_norm >= prev_res_norm:
                print('Iterative solver NOT CONVERGING!!! STOPPED!!!')
                return
            else:
                # print(f'Converge rate: {res_norm / prev_res_norm}')
                pass
            if self.CheckResidual(res_norm):
                return



