import os
import sys
import logging
import datetime
import numpy as np
import csv
from time import time, ctime, localtime, strftime

try:
    import hashlib
    md5_constructor = hashlib.md5
except ImportError:
    import md5
    md5_constructor = md5.new

__version__ = '$Id: files.py 685 2012-03-29 04:22:32Z carthur $'

logger = logging.getLogger()


def flModulePath(level=1):
    """
    Get the path of the module <level> levels above this function

    Input: level - level in the stack of the module calling this function
           (default = 1, function calling flModulePath)
    Output: path, basename and extension of the file containing the module
    Example: path, base, ext = flModulePath( )
             Calling flModulePath() from "/foo/bar/baz.py" produces the result
             "/foo/bar", "baz", ".py"
    """
    filename = os.path.realpath(sys._getframe(level).f_code.co_filename)
    path, fname = os.path.split(filename)
    base, ext = os.path.splitext(fname)
    path = path.replace(os.path.sep, '/')
    return path, base, ext


def flModuleName(level=1):
    """
    Get the name of the module <level> levels above this function

    Input: level - level in the stack of the module calling this function
           (default = 1, function calling flModuleName)
    Output: module name (as str)
    Example: mymodule = flModuleName( )
    """
    package = sys._getframe(level).f_code.co_name
    return package


def flProgramVersion(level=None):
    """
    Function to return the __version__ string from the parent
    program, where it is defined.
    If it is not defined, return an empty string.

    Input: level - level in the stack of the main script
           (default = maximum level in the stack)
    Output: version string (defined as the __version__ global variable)

    Example: my_program_version = flProgramVersion( )
    """
    if not level:
        import inspect
        level = len(inspect.stack()) - 1
    f = sys._getframe(level)
    if f.f_globals.has_key('__version__'):
        return f.f_globals['__version__']
    else:
        return ''


def flLoadFile(filename, comments='%', delimiter=',', skiprows=0):
    return np.genfromtxt(filename, comments=comments, delimiter=delimiter,
            skiprows=skiprows)


def flSaveFile(filename, data, header='', delimiter=' ', fmt='%.18e'):

    directory, fname = os.path.split(filename)
    if not os.path.isdir(directory):
        os.makedirs(directory)

    try:
        np.savetxt(filename, data, header=header, delimiter=delimiter, fmt=fmt)
    except TypeError:
        np.savetxt(filename, data, delimiter=delimiter, fmt=fmt)


def flGetStat(filename, CHUNK=2 ** 16):
    """
    Get basic statistics of filename - namely directory, name (excluding
    base path), md5sum and the last modified date. Useful for checking
    if a file has previously been processed.

    Input: filename, chunk size (for md5sum calculation)
    Output: path, name, md5sum, modification date
    Example: dir, name, md5sum, moddate = flGetStat(filename)
    """
    try:
        fh = open(filename)
        fh.close()
    except:
        logger.exception("Cannot open %s" % (filename))
        raise IOError("Cannot open %s" % (filename))

    try:
        directory, fname = os.path.split(filename)
    except:
        logger.exception('Input file is not a string')
        raise TypeError('Input file is not a string')

    try:
        si = os.stat(filename)
    except IOError:
        logger.exception('Input file is not a valid file: %s' % (filename))
        raise IOError('Input file is not a valid file: %s' % (filename))

    moddate = ctime(si.st_mtime)
    m = md5_constructor()
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

    Input: extension (default=.ini), prefix (default is empty)
    Output: Full path of calling function/module, with the source file's
    extension replaced with extension, and optionally prefix inserted
    after the last path separator
    Example: configFile = flConfigFile('.ini')
             Calling flConfigFile from /foo/bar/baz.py should
             return /foo/bar/baz.ini
    """
    if not level:
        import inspect
        level = len(inspect.stack())

    path, base, ext = flModulePath(level)
    config_file = os.path.join(path, prefix + base + extension)
    config_file = config_file.replace(os.path.sep, '/')
    return config_file


def flStartLog(logFile, logLevel, verbose=False, datestamp=False, newlog=True):
    """
    Start logging to logFile all messages of logLevel and higher.
    Setting verbose=True will report all messages to STDOUT as well

    Input: logFile - full path to log file
           logLevel - string specifiying one of the standard Python logging
               levels ('NOTSET','DEBUG','INFO','WARNING','ERROR','CRITICAL')
           verbose - boolean: True will echo all logging calls to STDOUT
           datestamp - boolean: True will include a timestamp of the creation
                       time in the filename
           newlog - boolean: True will create a new log file each time this
                    function is called. False will append to the existing file.
    Output: None
    Example: flStartLog('/home/user/log/app.log','INFO',verbose=True)
    """
    if datestamp:
        b, e = os.path.splitext(logFile)
        curdate = datetime.datetime.now()
        curdatestr = curdate.strftime('%Y%m%d%H%M')
        # The lstrip on the extension is required as splitext leaves it on.
        logFile = "%s.%s.%s" % (b, curdatestr, e.lstrip('.'))

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
    logger = logging.getLogger()

    if len(logger.handlers) < 2:
        # Assume that the second handler is a StreamHandler for verbose
        # logging. This ensures we do not create multiple StreamHandler
        # instances that will *each* print to STDOUT
        if verbose and sys.stdout.isatty():
            # If set to true, all logging calls will also be printed to the
            # console (i.e. STDOUT)
            console = logging.StreamHandler()
            console.setLevel(getattr(logging, logLevel))
            formatter = logging.Formatter('%(asctime)s: %(message)s',
                                          '%H:%M:%S', )
            console.setFormatter(formatter)
            logger.addHandler(console)

    logger.info('Started log file %s (detail level %s)' % (logFile, logLevel))
    logger.info('Running %s (pid %d)' % (sys.argv[0], os.getpid()))
    logger.info('Version %s' % (flProgramVersion()))
    return logger


def flLogFatalError(tblines):
    """
    Log the error messages normally reported in a traceback
    so that all error messages can be caught.
    the input 'tblines' is created by calling
    traceback.format_exc().splitlines()
    """
    for line in tblines:
        logger.critical(line.lstrip())
    sys.exit()


def flModDate(filename, dateformat='%Y-%m-%d %H:%M:%S'):
    """
    Return the update date of the input file

    Input: filename - file name (full path)
           dateformat - (optional) format string for the date
    Output: File modification date/time as a string
    Example: modDate = flModDate( 'C:/foo/bar.csv' , dateformat='%Y-%m-%dT%H:%M:%S' )
    """
    try:
        si = os.stat(filename)
    except IOError:
        logger.exception('Input file is not a valid file: %s' % (filename))
        raise IOError('Input file is not a valid file: %s' % (filename))
    moddate = localtime(si.st_mtime)

    return strftime(dateformat, moddate)


def flSize(filename):
    """
    Return the size of the input file in bytes

    Input: filename - file name (full path)
    Output: file size in bytes
    Example: file_size = flSize( 'C:/foo/bar.csv' )
    """
    try:
        si = os.stat(filename)
    except WindowsError:
        logger.exception('Input file is not a valid file: %s' % (filename))
        raise IOError('Input file is not a valid file: %s' % (filename))
    except IOError:
        logger.exception('Input file is not a valid file: %s' % (filename))
        raise IOError('Input file is not a valid file: %s' % (filename))
    else:
        size = si.st_size

    return size
