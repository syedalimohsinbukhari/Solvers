# Newton Interpolators

Named after their inventor, [Isaac Newton](https://en.wikipedia.org/wiki/Isaac_Newton), these methods make use of
the [divided difference method](https://en.wikipedia.org/wiki/Divided_differences) to calculate the coefficients of the
polynomial representing the given set of data points. The resulting polynomial is than used to calculate the value at
any given point. In case of this library, the value of $y(x_{\text{new}})$ is directly calculated.

The library includes implementation of,

* Forward (`FWD`) interpolation,
* Backward (`BKW`) interpolation,
* Divided difference (`DVd`) interpolation.

All three methods make use of a difference table, however, the method of difference table calculation is different for `DVd` than the other two methods.

## Calculation of difference table

For `FWD` and `BKW` methods, the difference table can be written as,

|  $x$  |  $f(x)$  |   $\Delta f_0$    |         $\Delta^2 f_0$          |
| :---: | :------: | :---------------: | :-----------------------------: |
| $x_0$ | $f(x_0)$ |                   |                                 |
|       |          | $f(x_1) - f(x_0)$ |                                 |
| $x_1$ | $f(x_1)$ |                   | $\Delta f(x_1) - \Delta f(x_0)$ |
|       |          | $f(x_2) - f(x_1)$ |                                 |
| $x_2$ | $f(x_2)$ |                   |                                 |

> [!NOTE]
> The notation for `FWD` and `BKW` differs, the `FWD` method uses $\Delta^{n-1} f_0$ whereas the `BKW` method uses $\nabla^{n-1} f_0$, but they represent the term in the table.

However, for `DVd` the method of calculating the difference table is a bit different

|  $x$  |  $f(x)$  |                $\Delta^\prime f_0$                |                          $\Delta^{\prime\prime} f_0$                          |
| :---: | :------: | :-----------------------------------------------: | :---------------------------------------------------------------------------: |
| $x_0$ | $f(x_0)$ |                                                   |                                                                               |
|       |          | $\displaystyle \frac{f(x_1) - f(x_0)}{x_1 - x_0}$ |                                                                               |
| $x_1$ | $f(x_1)$ |                                                   | $\displaystyle \frac{\Delta^\prime f(x_1) - \Delta^\prime f(x_0)}{x_2 - x_0}$ |
|       |          | $\displaystyle \frac{f(x_2) - f(x_1)}{x_2 - x_1}$ |                                                                               |
| $x_2$ | $f(x_2)$ |                                                   |                                                                               |

## Forward interpolation

For `FWD` interpolation, the following relationship is used,

$$f(x) = f(x_0) + p\Delta f_0 + \frac{p(p-1)}{2!}\Delta^2f_0 + \frac{p(p-1)(p-2)}{3!}\Delta^3f_0 + \cdots$$

## Backward interpolation

For `BKW` interpolation, the following relationship is used,

$$f(x) = f(x_0) + p\nabla f_0 + \frac{p(p+1)}{2!}\nabla f_0 + \frac{p(p+1)(p+2)}{3!}\nabla f_0 + \cdot$$

## Divided difference interpolation

For `DVd` interpolation, the following relationship is used,

$$f(x) = f(x_0) + (x-x_0)\Delta^\prime f_0 + (x-x_0)(x-x_1)\Delta^{\prime\prime}f_0 + \cdots$$












