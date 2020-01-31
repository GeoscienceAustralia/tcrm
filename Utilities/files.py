import os
import sys
import logging

import datetime
import numpy as np
from time import ctime, localtime, strftime

import hashlib

LOGGER = logging.getLogger()

if not getattr(__builtins__, "WindowsError", None):
    class WindowsError(OSError):
        pass

def flModulePath(level=1):
    """
    Get the path of the module <level> levels above this function

    :param int level: level in the stack of the module calling this function
                      (default = 1, function calling ``flModulePath``)

    :returns: path, basename and extension of the file containing the module

    Example: path, base, ext = flModulePath( )
    Calling flModulePath() from "/foo/bar/baz.py" produces the result
    "/foo/bar", "baz", ".py"
    """
    filename = os.path.realpath(sys._getframe(level).f_code.co_filename)
    path, fname = os.path.split(filename)
    path.replace(os.path.sep, '/')
    base, ext = os.path.splitext(fname)
    return path, base, ext


def flModuleName(level=1):
    """
    Get the name of the module <level> levels above this function

    :param int level: Level in the stack of the module calling this function
                      (default = 1, function calling ``flModuleName``)
    :returns: Module name.
    :rtype: str

    Example: mymodule = flModuleName( )
    Calling flModuleName() from "/foo/bar/baz.py" returns "baz"

    """
    package = sys._getframe(level).f_code.co_name
    return package


def flProgramVersion(level=None):
    """
    Return the __version__ string from the top-level program, where defined.

    If it is not defined, return an empty string.

    :param int level: level in the stack of the main script
                      (default = maximum level in the stack)
    :returns: version string (defined as the ``__version__`` global variable)

    """
    if not level:
        import inspect
        level = len(inspect.stack()) - 1
    f = sys._getframe(level)
    if '__version__' in f.f_globals:
        return f.f_globals['__version__']
    else:
        return ''


def flLoadFile(filename, comments='%', delimiter=',', skiprows=0):
    """
    Load a delimited text file -- uses :func:`numpy.genfromtxt`

    :param filename: File, filename, or generator to read
    :type  filename: file or str
    :param comments: (default '%') indicator
    :type  comments: str, optional
    :param delimiter: The string used to separate values.
    :type  delimiter: str, int or sequence, optional

    """
    return np.genfromtxt(filename, comments=comments,
                         delimiter=delimiter,
                         skip_header=skiprows)


def flSaveFile(filename, data, header='', delimiter=',', fmt='%.18e'):
    """
    Save data to a file.

    Does some basic checks to ensure the path exists before attempting
    to write the file. Uses :class:`numpy.savetxt` to save the data.

    :param str filename: Path to the destination file.
    :param data: Array data to be written to file.
    :param str header: Column headers (optional).
    :param str delimiter: Field delimiter (default ',').
    :param str fmt: Format statement for writing the data.

    """

    directory, fname = os.path.split(filename)
    if not os.path.isdir(directory):
        os.makedirs(directory)

    try:
        np.savetxt(filename, data, header=header, delimiter=delimiter, fmt=fmt,
                   comments='%')
    except TypeError:
        np.savetxt(filename, data, delimiter=delimiter, fmt=fmt, comments='%')


def flGetStat(filename, CHUNK=2 ** 16):
    """
    Get basic statistics of filename - namely directory, name (excluding
    base path), md5sum and the last modified date. Useful for checking
    if a file has previously been processed.

    :param str filename: Filename to check.
    :param int CHUNK: (optional) chunk size (for md5sum calculation).

    :returns: path, name, md5sum, modification date for the file.
    :raises TypeError: if the input file is not a string.
    :raises IOError: if the file is not a valid file, or if the file
                     cannot be opened.

    Example: dir, name, md5sum, moddate = flGetStat(filename)

    """
    try:
        fh = open(filename)
        fh.close()
    except:
        LOGGER.exception("Cannot open %s" % (filename))
        raise IOError("Cannot open %s" % (filename))

    try:
        directory, fname = os.path.split(filename)
    except:
        LOGGER.exception('Input file is not a string')
        raise TypeError('Input file is not a string')

    try:
        si = os.stat(filename)
    except IOError:
        LOGGER.exception('Input file is not a valid file: %s' % (filename))
        raise IOError('Input file is not a valid file: %s' % (filename))

    moddate = ctime(si.st_mtime)
    m = hashlib.md5()
    f = open(filename, 'rb')

    while 1:
        chunk = f.read(CHUNK)
        if not chunk:
            break
        m.update(chunk)
    md5sum = m.hexdigest()

    return directory, fname, md5sum, moddate


def flConfigFile(extension='.ini', prefix='', level=None):
    """
    Build a configuration filename (default extension .ini) based on the
    name and path of the function/module calling this function. Can also
    be useful for setting log file names automatically.
    If prefix is passed, this is preprended to the filename.

    :param str extension: file extension to use (default '.ini'). The
                          period ('.') must be included.
    :param str prefix: Optional prefix to the filename (default '').
    :param level: Optional level in the stack of the main script
                  (default = maximum level in the stack).
    :returns: Full path of calling function/module, with the source file's
              extension replaced with extension, and optionally prefix
              inserted after the last path separator.

    Example: configFile = flConfigFile('.ini')
    Calling flConfigFile from /foo/bar/baz.py should return /foo/bar/baz.ini

    """

    if not level:
        import inspect
        level = len(inspect.stack())

    path, base, ext = flModulePath(level)
    config_file = os.path.join(path, prefix + base + extension)
    return config_file


def flStartLog(logFile, logLevel, verbose=False, datestamp=False, newlog=True):
    """
    Start logging to logFile all messages of logLevel and higher.
    Setting ``verbose=True`` will report all messages to STDOUT as well.

    :param str logFile: Full path to log file.
    :param str logLevel: String specifiying one of the standard Python logging
                         levels ('NOTSET','DEBUG','INFO','WARNING','ERROR',
                         'CRITICAL')
    :param boolean verbose: ``True`` will echo all logging calls to STDOUT
    :param boolean datestamp: ``True`` will include a timestamp of the creation
                              time in the filename.
    :param boolean newlog: ``True`` will create a new log file each time this
                           function is called. ``False`` will append to the
                           existing file.

    :returns: :class:`logging.logger` object.

    Example: flStartLog('/home/user/log/app.log', 'INFO', verbose=True)
    """
    if datestamp:
        base, ext = os.path.splitext(logFile)
        curdate = datetime.datetime.now()
        curdatestr = curdate.strftime('%Y%m%d%H%M')
        # The lstrip on the extension is required as splitext leaves it on.
        logFile = "%s.%s.%s" % (base, curdatestr, ext.lstrip('.'))

    logDir = os.path.dirname(os.path.realpath(logFile))
    if not os.path.isdir(logDir):
        try:
            os.makedirs(logDir)
        except OSError:
            # Unable to create the directory, so stick it in the
            # current working directory:
            path, fname = os.path.split(logFile)
            logFile = os.path.join(os.getcwd(), fname)

    if newlog:
        mode = 'w'
    else:
        mode = 'a'

    logging.basicConfig(level=getattr(logging, logLevel),
                        format='%(asctime)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename=logFile,
                        filemode=mode)
    LOGGER = logging.getLogger()

    if len(LOGGER.handlers) < 2:
        # Assume that the second handler is a StreamHandler for verbose
        # logging. This ensures we do not create multiple StreamHandler
        # instances that will *each* print to STDOUT
        if verbose:
            console = logging.StreamHandler(sys.stdout)
            console.setLevel(getattr(logging, logLevel))
            formatter = logging.Formatter('%(asctime)s: %(message)s',
                                          '%H:%M:%S', )
            console.setFormatter(formatter)
            LOGGER.addHandler(console)

    LOGGER.info('Started log file %s (detail level %s)' % (logFile, logLevel))
    LOGGER.info('Running %s (pid %d)' % (sys.argv[0], os.getpid()))
    return LOGGER


def flLogFatalError(tblines):
    """
    Log the error messages normally reported in a traceback so that
    all error messages can be caught, then exit. The input 'tblines'
    is created by calling ``traceback.format_exc().splitlines()``.

    :param list tblines: List of lines from the traceback.

    """
    for line in tblines:
        LOGGER.critical(line.lstrip())
    sys.exit()


def flModDate(filename, dateformat='%Y-%m-%d %H:%M:%S'):
    """
    Return the last modified date of the input file

    :param str filename: file name (full path).
    :param str dateformat: Format string for the date (default
                           '%Y-%m-%d %H:%M:%S')
    :returns: File modification date/time as a string
    :rtype: str

    Example: modDate = flModDate( 'C:/foo/bar.csv' , dateformat='%Y-%m-%dT%H:%M:%S' )
    """
    try:
        si = os.stat(filename)
    except (IOError, WindowsError):
        LOGGER.exception('Input file is not a valid file: %s' % (filename))
        raise IOError('Input file is not a valid file: %s' % (filename))
    moddate = localtime(si.st_mtime)

    return strftime(dateformat, moddate)


def flSize(filename):
    """
    Return the size of the input file in bytes

    :param str filename: Full path to the file.
    :returns: File size in bytes.
    :rtype: int

    Example: file_size = flSize( 'C:/foo/bar.csv' )
    """
    try:
        si = os.stat(filename)
    except (IOError, WindowsError):
        LOGGER.exception('Input file is not a valid file: %s' % (filename))
        raise OSError('Input file is not a valid file: %s' % (filename))
    else:
        size = si.st_size

    return size
