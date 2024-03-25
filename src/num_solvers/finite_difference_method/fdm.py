"""Created on Jan 11 13:59:35 2024"""

import numpy as np

from src.num_solvers.__backend.fdm_ import OneDimensionalFDM, OneDimensionalPDESolver

L = 2
x_range, dx, dt, k = [0, L], 0.1, 0.01, 9
x_r = np.linspace(0, L, int(L / dx))

oneD_FDM = OneDimensionalFDM(x_range, dx, dt)

d = OneDimensionalPDESolver(lambda x: 2 * x**2 - 4 * x,
                            [0, 0],
                            [oneD_FDM.d1_backward(), -k * oneD_FDM.d2_central()],
                            ic_values=x_r)

d1 = d.solve(3)

print(d1[0])
print(d1[1])

# def exact_solution(n, i_dx, i_dt, k):
#     sum_ = 0
#     for i in range(n):
#         if i % 2 == 1:
#             f1_ = 400 / (pi**3 * i**3)
#             f1_ *= sin(i * pi * i_dx) * exp(-i**2 * pi**2 * k * i_dt)
#             sum_ += f1_

#     return sum_


# plt.plot(x_, [i - j for i, j in zip(d1[0], [exact_solution(20, i, 1, k) for i in x_])],
#          label='numerical - original n=20 difference')
# plt.plot(x_, [i - j for i, j in zip(d1[0], [exact_solution(200, i, 1, k) for i in x_])],
#          label='numerical - original n=200 difference')
# plt.plot(x_, [i - j for i, j in zip(d1[0], [exact_solution(2000, i, 1, k) for i in x_])],
#          label='numerical - original n=2000 difference')

# plt.plot(x_, [i - j for i, j in zip(d1[1], [exact_solution(20, i, 2, k) for i in x_])],
#          label='numerical - original n=20 difference')
# plt.plot(x_, [i - j for i, j in zip(d1[1], [exact_solution(200, i, 2, k) for i in x_])],
#          label='numerical - original n=200 difference')
# plt.plot(x_, [i - j for i, j in zip(d1[1], [exact_solution(2000, i, 2, k) for i in x_])],
#          label='numerical - original n=2000 difference')

# plt.legend(loc='best')
# plt.tight_layout()
# plt.show()
