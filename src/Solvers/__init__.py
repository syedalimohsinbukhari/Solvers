"""Created on Dec 20 14:10:28 2023"""

from typing import Callable, Optional, Union

from numpy.typing import NDArray

TOLERANCE = 1e-10

DOC_STYLE = 'numpy_napoleon_with_merge'

IFloat = Union[float, int]

FList = list[IFloat]
LList = list[FList]
LLList = list[LList]
OptList = Optional[FList]

FListOrLList = Union[FList, LList]
IFloatOrFList = Union[IFloat, FList]

OptFunc = Optional[Callable]
LFunc = list[Callable]

NdArray2 = tuple[NDArray, NDArray]
