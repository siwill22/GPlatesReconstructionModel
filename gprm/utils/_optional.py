"""Helpers for optional dependencies.

Several gprm functions need a package that is not part of the core install, usually because
it is awkward to install (pygmt needs the GMT command-line library, which pip cannot provide)
and only a minority of the library needs it. Those functions import it at the point of use
and call :func:`require` so that a missing dependency says which extra to install rather than
just naming the module.
"""

# Which extra installs what, for the error message.
_INSTALL_HINT = {
    'pygmt': ('viz', 'pygmt also needs the GMT command-line library, which pip cannot install '
                     'for you: `conda install -c conda-forge gmt`, Homebrew\'s `gmt`, or your '
                     'system package manager'),
    'stripy': ('spatial', None),
    'litho1pt0': ('geophysics', None),
    'pmagpy': ('geophysics', None),
    'pyshtools': ('geophysics', None),
}


def require(module_name, used_for=None):
    """Import an optional dependency, or raise an ImportError that says how to get it.

    :param module_name: Name of the module to import.
    :param used_for: Short description of what needs it, used in the message.
    :returns: The imported module.
    """
    try:
        return __import__(module_name)
    except ImportError as error:
        extra, note = _INSTALL_HINT.get(module_name, (None, None))

        message = "{:s} is required{:s}, but is not installed.".format(
            module_name, ' for {:s}'.format(used_for) if used_for else '')
        if extra:
            message += "\nInstall it with:  pip install 'gprm[{:s}]'".format(extra)
        if note:
            message += "\nNote: {:s}.".format(note)

        raise ImportError(message) from error
