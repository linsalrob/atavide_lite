"""
An Error Class so I can write my own errors
"""
class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass

class ColorNotFoundError(Error):
    """
    Exception raised for a color not being found.

    :param message: explanation of the error
    """

    def __init__(self, message):
        self.message = message
