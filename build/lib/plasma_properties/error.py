class Error(Exception):
    """Base class for other exceptions"""
    pass

class DimensionMismatchError(Error):
    """Raised when arrays don't match size"""
    pass

class ValueError(Error):
    """Raised when value is incorrect"""
    pass

class TypeError(Error):
    """Raised when type is incorrect"""
    pass