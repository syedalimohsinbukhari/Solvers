# Cholesky decomposition

Discovered by [André-Louis Cholesky](https://en.wikipedia.org/wiki/Andr%C3%A9-Louis_Cholesky), this method is applicable for [Hermitian](https://en.wikipedia.org/wiki/Hermitian_matrix), symmetric and [positive-definite](https://en.wikipedia.org/wiki/Positive-definite_matrix) matrices.

For each applicable matrix, the Cholesky decomposition reduces the given matrix into a lower-triangular matrix and its conjugate transpose matrix.

Mathematically, it can be written as,

$$\textbf{A} = \textbf{L}\textbf{L}^*$$

For a non-complex matrix, the above statement can also be written as,

$$\textbf{A} = \textbf{L}\textbf{L}^T$$

> [!NOTE]
> The Cholesky decomposition can only be performed only on positive-definite matrices.

## Method

For a matrix,

$$
A = 
\begin{bmatrix}
    a_{11} & b_{12} & c_{13}\\
    b_{12} & b_{22} & c_{23}\\
    c_{13} & c_{23} & c_{33}
\end{bmatrix}
$$

The following formulae can be used to construct the lower tri-angular matrix,

$$
\ell_{ki} = \frac{1}{\ell_{ii}}\left(a_{ki} - \sum_{j=1}^{i-1}\ell_{ij}\cdot\ell_{kj}\right)
$$

and

$$
\ell_{kk} = \sqrt{a_{kk} - \sum_{j=1}^{k-1}\ell_{kj}^2}
$$

which form the elements of matrix $L$ such that,

$$
L = 
\begin{bmatrix}
    \ell_{11} & 0 & 0\\
    \ell_{21} & \ell_{22} & 0\\
    \ell_{31} & \ell_{32} & \ell_{33}
\end{bmatrix}
$$

The other matrix can be obtained by simply transposing the matrix $L$.