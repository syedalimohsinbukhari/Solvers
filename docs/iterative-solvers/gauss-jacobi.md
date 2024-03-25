# Gauss Jacobi iterative method

In numerical linear algebra, the Jacobi method (a.k.a. the Jacobi iteration method) is an iterative algorithm for
determining the solutions of
a [strictly diagonally dominant](https://en.wikipedia.org/wiki/Diagonally_dominant_matrix) [system of linear equations](https://en.wikipedia.org/wiki/System_of_linear_equations).
Each diagonal element is solved for, and an approximate value is plugged in. The process is then iterated until it
converges.

## Method

The element-based formula for solution of system

$$
\begin{equation}
x_i^{k+1} = \frac{1}{a_{ii}}\left(b_i - \sum_{j\neq i}a_{ij}x_k^{(k)}\right)
\end{equation}
$$

where $i = 1, 2, 3, \ldots, n$. The procedure is generally continued until the changes made by an iteration are below
some tolerance, such as a sufficiently small residual. It should be noted that the Gauss-Jacobi method is very sensitive
and only works well with strictly diagonally dominant matrices, such that,

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
> See alo, *
*[Jacobi-Matrix method](https://github.com/syedalimohsinbukhari/Solvers/blob/md-docs/docs/matrix-decomposition/gauss-jacobi.md)
**.