# Lagrange Interpolator

Similar to the Newton interpolator, the Lagrange interpolator method is used to estimate the polynomial for the given data sets. Once the polynomial has been determined, the interpolation can be performed. The Lagrange interpolation formula can be written as,

$$L(x) = \displaystyle \sum_{j=0}^k f(x_j)\left(\prod_{0\leq m \leq k}^{m\neq j} \frac{x-x_m}{x_j-x_m}\right)$$


