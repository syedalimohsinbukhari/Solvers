# Spline interpolation

A form of interpolation which deploys piecewise functions between all given points for interpolation. Instead of fitting
a single, high-degree polynomial, the spline interpolation fits low-degree polynomials between all points.

The three types of spline interpolation implemented in this library are,

1. Linear spline
2. Quadratic spline
3. Natural cubic spline

## Linear spline

The linear spline method tries to fit the data using a set of line segments between two adjacent data points. This can
be written mathematically as,

$$ S_i(x) = \frac{x - x_i}{x_{i+1} - x_i}f(x_{i+1}) - \frac{x-x_{i+1}}{x_{i+1} - x_i}f(x_i) $$

Where $S_i(x) \in [x_i, x_{i+1}]$

## Quadratic spline

The quadratic spline method fits quadratic equations between adjacent data points. The method of making quadratic spline
comprises several steps,

### How many equations?

As stated earlier, the quadratic splines make use of quadratic equations for fitting the adjacent data points, this
means that for each spline constructed, we have three unknown and thus ideally should have three equations,

$$ a_1x_1^2 + b_1x_1+c_1=0$$

Consider a data set with four data points, $\langle x_1, y_1\rangle, \langle x_2, y_2\rangle, \langle x_3, y_3\rangle$
and $\langle x_4, y_4\rangle$. So for four data points, we required $3(N-1)$ equations.

Now consider the following table,

|   x   |   1   |   2   |   3   |   4   |
| :---: | :---: | :---: | :---: | :---: |
| y(x)  |   7   |  16   |   5   |   8   |

### Step 1:

Construct splines between each adjacent point,

$$a_1x_1^2 + b_1x_1 + c_1 = y(x_1)\quad;\quad 0 \leq x \leq 1$$

$$a_2x_2^2 + b_2x_2 + c_2 = y(x_2)\quad;\quad 1 \leq x \leq 2$$

$$a_3x_3^2 + b_3x_3 + c_3 = y(x_3)\quad;\quad 2 \leq x \leq 3$$

Which gives us,

$$ a_1x_1^2 + b_1x_1 + c_1 = 7$$

$$ a_1x_2^2 + b_1x_2 + c_1 = 16$$

$$ a_2x_2^2 + b_2x_2 + c_2 = 16$$

$$ a_2x_3^2 + b_2x_3 + c_2 = 5$$

$$ a_3x_3^2 + b_3x_3 + c_3 = 5$$

$$ a_3x_4^2 + b_3x_4 + c_3 = 8$$

Now we have $6$ equations for $9$ unknown variables. More equations can be developed at mid-point $(1, 2)$, and $(2, 3)$
by taking the derivatives of the spline equations at those points.

$$ 2a_1x + b_1 -2a_2x - b_2 = 0$$

$$ 2a_2x + b_2 -2a_3x - b_3 = 0$$

Now we have $8$ equations, the final equation is taken as $a_{N} = 0$, or in our case, $a_4 = 0$. This is because we
assume that the line leaving the endpoint is a straight line, which zeros the quadratic component.

Now, having all $9$ equations, we can write them in form of a matrix,

$$
\begin{bmatrix}
1 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
4 & 2 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 4 & 2 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 9 & 3 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 9 & 3 & 1\\
0 & 0 & 0 & 0 & 0 & 0 & 16 & 4 & 1\\
4 & 1 & 0 & -4 & -1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 6 & 1 & 0 & -6 & -1 & 0\\
1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
a_1\\\ b_1\\\ c_1\\\ a_2 \\\ b_2 \\\ c_2 \\\ a_3 \\\ b_3 \\\ c_3
\end{bmatrix} =
\begin{bmatrix}
7 \\\ 16 \\\ 16 \\\ 5 \\\ 5 \\\ 8 \\\ 0 \\\ 0 \\\ 0
\end{bmatrix}
$$

This matrix system can be solved using any numerical technique, following

$$ A\mathcal{x} = B \rightarrow \mathcal{x} = A^{-1}B$$

or with any other matrix decomposition, e.g., QR decomposition,

$$ A\mathcal{x} = B \rightarrow QR\mathcal{x} = B \rightarrow \mathcal{x} = R^{-1}Q^{-1}B$$

The resulting solution matrix, for this problem is,

$$
\text{solution} =
\begin{bmatrix}
0 & 9 & -2.0 & -20.0 & 89.0 & -82.0 & 34.0 & -235.0 & 404.0
\end{bmatrix}^T
$$

## Natural-cubic spline

As the name suggests, the natural-cubic spline uses a cubic equation between adjacent points. The method of developing a
cubic spline in similar to quadratic spline. The main difference being the evaluation of cubic equation instead of
quadratic. To form extra equations, form the inside points, the natural-cubic spline requires second and third
derivative as well.