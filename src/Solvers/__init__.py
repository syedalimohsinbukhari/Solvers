"""Created on Dec 20 14:10:28 2023"""

from typing import Callable, List, Optional, Union

TOLERANCE = 1e-10

DOC_STYLE = 'numpy_napoleon_with_merge'

IFloat = Union[float, int]

FList = List[IFloat]
LList = List[FList]
LLList = List[LList]
OptList = Optional[FList]

FListOrLList = Union[FList, LList]
IFloatOrFList = Union[IFloat, FList]

OptFunc = Optional[Callable]
