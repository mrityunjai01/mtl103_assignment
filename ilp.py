import numpy as np
import pdb
# check if no zero div

def gomory(fn):
    with open(fn) as f:
        n, m = [int(x) for x in next(f).split()]
        b = np.array([int(x) for x in next(f).split()])
        c = np.array([int(x) for x in next(f).split()])
        a = []
        for i in range(m):
            a.append([int(x) for x in next(f).split()])

        a = np.array(a)

        assert(a.shape == (m, n))
        assert(b.shape == (m,))
        assert(c.shape == (n,))
        lp(a, b, c, m, n)

class Problem:
    def __init__ (self, A, b, c, m, n, slack = False, aux = False):
        self.A = A # m X n
        self.b = b # m
        self.c = c # n
        self.slack = slack
        self.m = m
        self.eps = 1e-15
        self.n = n
        self.aux = aux
        if (self.slack):
            self.A = np.c_[self.A, np.eye(m)]
            self.c = np.r_[self.c, np.zeros(m)]
            assert(self.A.shape == (self.m, self.n + self.m))
            assert(self.b.shape == (self.m, ))
            # initial basis is n, ..., n + m - 1
            self.basis = np.array(list(range(n, n+m)))
            self.cbt = np.zeros((1, m))
            if (self.aux):

                for i in range(m):
                    if (self.b[i] < 0):
                        self.b[i] = -self.b[i]
                        self.A[i, :] = -self.A[i, :]
                        self.A = np.c_[self.A, ]
                        self.
    def solve(self):
        self.b_ = self.b.astype('float64')
        iteration = 0
        self.r = self.c.T - self.cbt @ self.A
        self.r = self.r.flatten()
        # pdb.set_trace()
        n, m = self.n, self.m
        while True:
            print(f"iteration {iteration}")
            incoming = np.argmax(self.r)
            max_r = np.max(self.r)
            print(f"max_r {max_r}")
            if (max_r <= self.eps):
                break
            ratios = np.where(self.A[:, incoming] > 0, self.b_ / (self.A[:, incoming] + self.eps), np.Infinity)
            print(ratios)
            outgoing = np.argmin(ratios)


            temp_basis = self.basis.tolist()

            temp_basis.pop(outgoing)

            temp_basis.append(incoming)

            self.basis = np.array(temp_basis)

            assert(self.basis.shape == (self.m,))

            self.basis_mask = np.zeros(n + m,)

            np.put(self.basis_mask, self.basis, 1)

            self.cbt = self.c[self.basis]

            # outgoing row, incoming col
            for i in range(m):

                if (i==outgoing):
                    self.b_ /= self.A[outgoing, incoming]
                    self.A[outgoing, :] /= self.A[outgoing, incoming]
                    continue


                assert(self.A[outgoing, incoming] != 0)
                self.A[i, :] = self.A[i, :] - self.A[i, incoming] / self.A[outgoing, incoming] * self.A[outgoing, :]

                self.b_[i] = self.b_[i] - self.A[i, incoming] / self.A[outgoing, incoming] * self.b_[outgoing]


            # pdb.set_trace()
            self.r = self.r - self.r[incoming] / self.A[outgoing, incoming] * self.A[outgoing, :]

            iteration = iteration + 1

        solution = np.zeros(n)
        solution[self.basis] = 1
        solution = np.where(solution > 0, self.b_, 0)
        basis =
        return solution, basis


def lp(A, b, c, m, n):
    # solve lp using two phase
    # phase 1
    p = Problem(A, b, c, m, n, slack = True)
    ans = p.solve()
    # phase 2
    p = Problem(A, b, c, m, n, slack = True)
    ans = p.solve()


def dual_simp(A, x, b):
    pass


# actually solve
if __name__ == '__main__':
    gomory('inp.txt')
