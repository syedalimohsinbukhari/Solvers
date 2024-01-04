"""Created on Dec 20 14:10:28 2023"""

from typing import Callable, List, Optional, Union

TOLERANCE = 1e-10

DOC_STYLE = 'numpy_napoleon_with_merge'

FLOAT_OR_INT = Union[float, int]
LIST = List

F_LIST = LIST[FLOAT_OR_INT]
L_LIST = LIST[F_LIST]
L_L_LIST = List[L_LIST]
O_LIST = Optional[F_LIST]

FLIST_OR_LLIST = Union[F_LIST, L_LIST]
FLINT_OR_FLIST = Union[FLOAT_OR_INT, F_LIST]

O_CALLABLE = Optional[Callable]
