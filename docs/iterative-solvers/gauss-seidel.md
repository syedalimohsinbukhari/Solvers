# Gauss Seidel iterative method

In numerical linear algebra, the Gauss-Seidel is an iterative algorithm for determining the solutions of a [strictly diagonally dominant](https://en.wikipedia.org/wiki/Diagonally_dominant_matrix) [system of linear equations](https://en.wikipedia.org/wiki/System_of_linear_equations). Each diagonal element is solved for, and an approximate value is plugged in. The process is then iterated until it converges. 

## Method

The element-based formula for solution of system of equations is,

$$
x_i^{k+1} = \frac{1}{a_{ii}}\left(b_i - \sum_{j=1}^{i-1}a_{ij}x_j^{k+1} - \sum_{j=i+1}^na_{ij}x_j^k\right)
$$

where $i = 1, 2, \ldots, n.$. The procedure is generally continued until the changes made by an iteration are below some tolerance, such as a sufficiently small residual. It should be noted that the like Gauss-Jacobi method, the Gauss-Seidel mehtod is also very sensitive and only works well with strictly diagonally dominant matrices, such that,

$$
\left|a_{11}\right| \geq \left|a_{12}\right| + \left|a_{13}\right|
$$

$$
\left|a_{22}\right| \geq \left|a_{21}\right| + \left|a_{23}\right|
$$

$$
\left|a_{33}\right| \geq \left|a_{31}\right| + \left|a_{32}\right|
$$



> [!NOTE]
> See alo, **[Gauss-Seidel matrix method](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/matrix-decomposition/gauss-seidel.md)**.