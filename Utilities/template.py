"""
:mod:`template` -- replace strings with another value
=====================================================

.. module:: template
    :synopsis: Read in a text file, with formatted elements
               and replace with strings.

.. moduleauthor:: Craig Arthur <craig.arthur@ga.gov.au>
"""
import re

def replace(infile, outfile, replacements):
    """
    Replace all instances of the keys with values in infile and w
    write to outfile.

    In the input file, keywords to be replaced should be written
    as '{keyword}'.

    e.g. if infile has a line containing a formatted keyword::

        Input = {PATH}

    and ``replacements = dict('PATH': '/foo/baz')``, then the
    output file will be written as::

        Input = /foo/baz


    :param str infile: Path to an input file.
    :param str outfile: Destination file to be written.
    :param dict replacements: A set of key-value pairs that dictate
              the replacements to occur in the input file.

    """
    fi = open(infile, 'r')
    fo = open(outfile, 'w')
    for line in fi:
        newline = line
        for key, val in list(replacements.items()):
            newline = re.sub('{'+key+'}', val, newline)

        fo.write(newline)
    fi.close()
    fo.close()
    return True
