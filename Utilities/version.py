"""
:mod:`version` -- provide details of software version
=====================================================

.. module:: version
    :synopsis: Provide ways to determine the version of
               the currently executing code.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>

"""

import os
import subprocess
from .files import flModulePath, flModDate


UPDATE_MSG = """
----------------------------------------------------------
Your TCRM version is not up-to-date. The last 3 things that
have been fixed are:

{0}

To update your version of TCRM, try running:

git pull

----------------------------------------------------------
"""

def git(command):
    """
    Execute the given command with git

    :param str command: A valid git command.
    :return: Output from the given command.
    :rtype: str

    """
    with open(os.devnull) as devnull:
        return subprocess.check_output('git ' + command,
                                       shell=True,
                                       stderr=devnull)
def status():
    """
    Check status TCRM of code.

    :returns: A message containing a listing of recent
              changes to the model code.
    :rtype: str

    """

    msg = ''

    try:
        git('fetch')  # update status of remote
        behind = int(git('rev-list HEAD...origin/master --count'))
        recent = git('log --pretty=format:" - %s (%ad)" --date=relative ' +
                     'origin/master HEAD~3..HEAD')
        if behind != 0:
            msg = UPDATE_MSG.format(recent)
    except subprocess.CalledProcessError:
        pass

    return msg

def version():
    """
    Check version of TCRM code.

    :returns: Current git commit hash and the date/time
             of the commit.

    :rtype: str
    :raises: :mod:`subprocess.CalledProcessError` if the git
             command fails.

    .. note:: This requires ``git`` to be installed to execute.
    """

    vers = ''

    try:
        vers = git('log -1 --date=iso --pretty=format:"%H"')
    except subprocess.CalledProcessError:
        # Case for missing git:
        import inspect
        path, name, ext = flModulePath(len(inspect.stack()))
        fname = os.path.join(path, name+ext)
        fdate = flModDate(fname)
        vers = '{0} modified {1}'.format(fname, fdate)

    return vers
