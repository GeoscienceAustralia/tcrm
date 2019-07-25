"""
:mod:`pathLocator` -- determine directory one level above the Utilites folder
=============================================================================

.. module:: pathLocator
    :synopsis: A function to determine the directory one level above
               the Utilities folder. Designed to work even if the code has been
               compiled with py2exe.  This results in the modules being built-in
               to the interpreter, so the path of the executable is returned
               instead.

.. moduleauthor:: Nicholas Summons <nicholas.summons@ga.gov.au>

"""
import sys
import os

def is_frozen():
    """
    Determine if modules have been built into the interpreter, e.g. by
    py2exe.

    :return: `True` if the modules are frozen, `False` otherwise.
    """
    # Modules built-in to the interpreter, e.g., by py2exe
    return hasattr(sys, "frozen")

def getRootDirectory():
    """
    Return the name of the path one level above the directory of this current
    file.

    :return: Path name one level above this.
    :rtype: str
    """

    encoding = sys.getfilesystemencoding()
    if is_frozen():
        return os.path.dirname(str(sys.executable, encoding))
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
