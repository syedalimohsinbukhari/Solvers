# Multi-RK implementation

The RK methods can be extended to be implemented on any number of ODEs, but this requires solving the ODEs simultaneously,

Consider two ODEs,

$$ y' = F(x, y, z) $$

$$ z' = G(x, y, z) $$

given that, $y(x_0) = y_0$ and $z(x_0) = z_0$, the implementation of multi-RK4 is as follows,

$$ y_{n+1} = y_n + \frac{1}{6}\left(k_1 + 2k_2 + 2k_3 + k_4\right)$$

$$ z_{n+1} = z_n + \frac{1}{6}\left(\ell_1 + 2\ell_2 + 2\ell_3 + \ell_4\right)$$

where solving the terms of $k_n$ and $\ell_n$ follows the following pattern,

$$ k_1 = hF(x_n, y_n, z_n)$$

$$ \ell_1 = hG(x_n, y_n, z_n)$$

$$ k_2 = hF\left(x_n + \frac{h}{2}, y_n + \frac{k_1}{2}, z_n + \frac{\ell_1}{2}\right)$$

$$ \ell_2 = hG\left(x_n + \frac{h}{2}, y_n + \frac{k_1}{2}, z_n + \frac{\ell_1}{2}\right)$$

$$ k_3 = hF\left(x_n + \frac{h}{2}, y_n + \frac{k_2}{2}, z_n + \frac{\ell_2}{2}\right)$$

$$ \ell_3 = hG\left(x_n + \frac{h}{2}, y_n + \frac{k_2}{2}, z_n + \frac{\ell_2}{2}\right)$$

$$ k_4 = hF\left(x_n + h, y_n + k_3, z_n + \ell_3 \right)$$

$$ \ell_4 = hG\left(x_n + h, y_n + k_3, z_n + \ell_3 \right)$$

This sequence of solving the equations is essential as each step requires the value for each preceeding RK-step to be solved. This can be extended to N-ODEs.