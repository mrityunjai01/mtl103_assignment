import numpy as np
import pdb
eps = 1e-15
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
        # pdb.set_trace()
        assert(a.shape == (m, n))
        assert(b.shape == (m,))
        assert(c.shape == (n,))
        ans = overall_wrapper(a, b, c, m, n)
        print("hre", ans)
        return ans


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
        self.an = 0
        if (self.slack):
            self.A = np.c_[self.A, np.eye(m)]
            self.c = np.r_[self.c, np.zeros(m)]
            assert(self.A.shape == (self.m, self.n + self.m))
            assert(self.b.shape == (self.m, ))
            # initial basis is n, ..., n + m - 1
            self.basis = list(range(n, n+m))
            self.cbt = np.zeros((1, m))
            if (self.aux):

                for i in range(m):
                    if (self.b[i] < 0):
                        self.b[i] = -self.b[i]
                        self.A[i, :] = -self.A[i, :]
                        self.A = np.c_[self.A, np.zeros(m)]
                        self.A[i, -1] = 1
                        self.basis.remove(n + i)
                        self.basis.append(n + m + self.an - 1)
                        self.an += 1
                self.c = np.r_[np.zeros(n + m), -np.ones(self.an)]

                # print(f"added aux, {self.A}, {self.b}, {self.c}")
    def solve(self, verbose = False):
        self.b_ = self.b.astype('float64')
        iteration = 0
        self.r = self.c.T - self.cbt @ self.A
        self.r = self.r.flatten()
        if (verbose):
            print(f"starting {self.A}, {self.b_}")
        # pdb.set_trace()
        n, m = self.n, self.m
        while True:
            # print(f"iteration {iteration}")
            incoming = np.argmax(self.r)
            max_r = np.max(self.r)
            # print(f"max_r {max_r}")
            if (max_r <= self.eps):
                break
            ratios = np.where(self.A[:, incoming] > 0, self.b_ / (self.A[:, incoming] + self.eps), np.Infinity)
            # print(ratios)
            outgoing = np.argmin(ratios)
            self.basis.pop(outgoing)
            self.basis.append(incoming)
            assert(len(self.basis) == self.m)
            self.basis_mask = np.zeros(n + m + self.an,)
            np.put(self.basis_mask, self.basis, 1)
            self.cbt = self.c[self.basis]

            # outgoing row, incoming col
            for i in range(m):

                if (i==outgoing):
                    self.b_[i] /= self.A[outgoing, incoming]
                    self.A[outgoing, :] /= self.A[outgoing, incoming]
                    continue


                assert(self.A[outgoing, incoming] != 0)

                self.b_[i] = self.b_[i] - self.A[i, incoming] / self.A[outgoing, incoming] * self.b_[outgoing]

                self.A[i, :] = self.A[i, :] - self.A[i, incoming] / self.A[outgoing, incoming] * self.A[outgoing, :]




            # pdb.set_trace()
            self.r = self.r - self.r[incoming] / self.A[outgoing, incoming] * self.A[outgoing, :]
            iteration = iteration + 1
            if (verbose):
                print(f"at iteration {iteration} {self.A}, {self.b_}")

        solution = np.zeros(n + m + self.an)
        solution[self.basis] = 1
        # solution = np.where(solution > 0, self.b_, 0)
        cost = np.dot(self.c, solution)
        return self.A, self.b_, self.c, self.r, cost


def lp(A, b, c, m, n):
    # solve lp using two phase
    # phase 1
    p = Problem(A, b, c, m, n, slack = True, aux = True)
    ans = p.solve()
    # print(f"got {ans}")

    cost_for_aux = ans[4]
    # if (cost_for_aux <= eps):
    #     print(f"aux okay")
    # else:
    #     print(f"aux not okay")
    #     return
    # phase 2
    p = Problem(A, b, c, m, n, slack = True)
    ans = p.solve()
    return ans

def frac_part(x):
    return np.where(x >= 0, np.modf(x)[0], 1 + np.modf(x)[0])
def is_not_integer(x):
    frac = frac_part(x)
    if (frac < 1 - eps and frac > eps):
        return True
    return False

def dual_simp(A, b, c, r, m, n):
    # print(f"starting dual simp")
    added_vars = 0
    b_list = b.tolist()
    for i, val in enumerate(b_list):
        if (val < -eps):

            # pdb.set_trace()
            incoming = np.argmin(np.where((A[i] != 0) & (r != 0), np.abs(r/A[i]), np.Infinity))
            # print(f"selected incoming {incoming}")
            # incoming col, ith row
            for j in range(A.shape[0]):
                if (j != i):
                    b[j] = b[j] - A[j, incoming] / A[i, incoming] * b[i]
                    A[j, :] = A[j, :] - A[j, incoming] / A[i, incoming] * A[i, :]
                else:
                    b[j] = b[j] / A[i, incoming]
                    A[j, :] = A[j, :] / A[i, incoming]

            r = r - r[incoming] / A[i, incoming] * A[i, :]
    return A, b, c, r



def overall_wrapper(A, b, c, init_m, init_n):
    A_ = A.copy()
    b_ = b.copy()
    c_ = c.copy()
    A_, b_, c_, r, cost = lp(A_, b_, c_, A.shape[0], A.shape[1])
    # print(f"right aft lpss to make A_ {A_} b {b_} r {r}")

    iteration = 0
    max_iterations = 3
    added_vars = 0
    while (iteration < max_iterations):
        b_list = b_.tolist()
        found = True
        non_integer_row = -1
        for i, val in enumerate(b_list):
            if (is_not_integer(val)):
                found = False
                non_integer_row = i
        if (found):

            break
        else:
            m, n = A_.shape
            i = non_integer_row
            b_ = np.r_[b_, -frac_part(b_[i])]
            A_ = np.c_[A_, np.zeros(m)]
            r = np.r_[r, np.zeros(1)]
            A_ = np.r_[A_, [-frac_part(A_[i, :])]]
            A_[-1, -1] = 1
            # print(f"added new vars to make A_ {A_} b {b_} r {r}")
            added_vars += 1
            if (A_.shape != (m+1, n+1)):
                pdb.set_trace()
            assert(A_.shape == (m+1, n+1))
            A_, b_, c_, r = dual_simp(A_, b_, c_, r, A.shape[0], A.shape[1],)
            # print(f"after dual sis to make A_ {A_} b {b_} r {r}")

        iteration += 1
    final_answer = []
    # pdb.set_trace()
    for i in range(init_n):
        s = A_[:, i].sum()
        if (np.abs(s - 1) > eps):
            final_answer.append(0)
        else:
            single_one = False
            position_one = -1
            for j in range(A_.shape[0]):
                if (np.abs(A_[j, i] - 1) < eps):
                    single_one = True
                    position_one = j
            if (single_one):
                final_answer.append(b_[position_one])
            else:
                final_answer.append(0)
    return final_answer








# actually solve
if __name__ == '__main__':
    gomory('inp.txt')
