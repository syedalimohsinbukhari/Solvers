"""Created on Dec 20 22:59:44 2023"""

__all__ = ['newton_raphson_solver', 'trapezoid_rule', 'simpson_13', 'simpson_38', 'boole', 'weddle_3h10',
           'weddle_41h140']

from .NewtonCotesIntegrators import (boole, newton_raphson_solver, simpson_13, simpson_38, trapezoid_rule, weddle_3h10,
                                     weddle_41h140)
