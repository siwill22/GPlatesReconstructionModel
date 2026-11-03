"""Shared download helper for the dataset loaders.

Every loader in this subpackage ultimately calls :func:`pooch.retrieve`. When a download
fails, pooch raises whatever ``requests`` or its own hash check raised, which surfaces to
the user as a bare traceback with no indication of *which* dataset failed, where its cache
lives, or whether the problem is theirs to fix.

This module wraps that one call. The loaders keep calling ``_retrieve(...)`` with pooch's
signature unchanged; they simply import it from here instead of from ``pooch``, so there is
one place that knows how to explain a failure.

MIT License

Copyright (c) 2017-2021 Simon Williams
"""

import inspect as _inspect
import os as _os
import shutil as _shutil
import zipfile as _zipfile
from pathlib import Path as _Path

from pooch import Unzip as _Unzip
from pooch import retrieve as _pooch_retrieve


class DatasetFetchError(RuntimeError):
    """A gprm dataset could not be downloaded or verified.

    Subclasses :class:`RuntimeError`, so code that already catches broad exceptions is
    unaffected. The original exception is always chained as ``__cause__``.
    """


# Frames in these modules are gprm's own download plumbing, not the loader the user called.
_INTERNAL = frozenset(['_fetch.py', '_remote_zip.py'])


def _calling_loader():
    """Name of the nearest frame outside the download plumbing, i.e. the loader the user called."""
    for frame in _inspect.stack()[1:]:
        if _Path(frame.filename).name not in _INTERNAL:
            return frame.function
    return None


def _diagnose(err, url):
    """Return an actionable hint for a failed retrieve, or None if we cannot classify it."""
    text = str(err)

    # pooch raises ValueError for a checksum mismatch, with both hashes in the message.
    if isinstance(err, ValueError) and 'hash' in text.lower():
        return (
            "The file downloaded but its checksum did not match the one recorded in gprm.\n"
            "Either the download was truncated, or the file has been changed upstream.\n"
            "Delete the cached copy and try again; if it still fails, the recorded hash\n"
            "in gprm is out of date and needs updating."
        )

    status = getattr(getattr(err, 'response', None), 'status_code', None)
    if status in (401, 403):
        return (
            "The server refused the request ({}). The host may require registration, or\n"
            "may be blocking automated downloads. Try fetching the URL in a browser and\n"
            "placing the file in the cache directory by hand."
        ).format(status)
    if status in (404, 410):
        return (
            "The server says this file no longer exists ({}). The URL recorded in gprm has\n"
            "gone stale — the dataset has moved or been withdrawn. This needs fixing in gprm;\n"
            "please report it."
        ).format(status)
    if status is not None and 500 <= status < 600:
        return (
            "The server returned an error ({}). This is usually temporary — try again later."
        ).format(status)

    # Connection-level failures. Checked by name so this module does not import requests.
    if type(err).__name__ in ('ConnectionError', 'ConnectTimeout', 'ReadTimeout', 'Timeout',
                              'SSLError', 'TooManyRedirects'):
        return (
            "Could not reach the server. Check your network connection, and whether the host\n"
            "above is up — several gprm datasets are hosted on personal academic domains."
        )

    return None


def fetch_error(err, url, path=None, extra=None):
    """Build the :class:`DatasetFetchError` for a failed download.

    Shared by :func:`retrieve` and by ``_remote_zip.retrieve_zip_member``, so that however a
    dataset is fetched, a failure reads the same way.

    :param err: the original exception, chained by the caller via ``raise ... from err``.
    :param url: the URL that was being fetched.
    :param path: the cache directory, if known.
    :param extra: optional extra lines (e.g. which member of an archive was wanted).
    """
    loader = _calling_loader()
    lines = [
        "Could not fetch {}.".format(
            "the dataset '{}'".format(loader) if loader else "a gprm dataset"),
        "",
        "  url   : {}".format(url),
    ]
    for line in (extra or []):
        lines.append("  {}".format(line))
    if path is not None:
        lines.append("  cache : {}".format(path))
    lines += ["", "{}: {}".format(type(err).__name__, err)]

    hint = _diagnose(err, url)
    if hint:
        lines += ["", hint]

    return DatasetFetchError("\n".join(lines))


def retrieve(url, known_hash, **kwargs):
    """:func:`pooch.retrieve` with a failure message that says what went wrong.

    Arguments are passed through unchanged. On success the return value is pooch's, so this
    is a drop-in replacement. On failure a :class:`DatasetFetchError` is raised naming the
    dataset, the URL and the cache directory, with the original exception chained.
    """
    try:
        return _pooch_retrieve(url, known_hash, **kwargs)
    except Exception as err:
        raise fetch_error(err, url, path=kwargs.get('path')) from err


class UnzipTidy(_Unzip):
    """:class:`pooch.Unzip` that leaves out unwanted archive members and renames others as it
    extracts them.

    pooch extracts an archive again whenever any of its members is missing from the extract
    folder. Loaders that deleted or renamed files after extracting (unwanted
    DeformingMeshPoints folders, ``.cache`` files, file names with spaces) therefore unpacked
    the whole archive again on every call, which is slow and gives every file a new
    modification time, so anything cached downstream on those times was thrown away. Doing the
    tidying here keeps the folder exactly as pooch expects it, so later calls extract nothing.

    :param extract_dir: as for :class:`pooch.Unzip`.
    :param skip: function of an archive member name (``'folder/file.gpml'``) that is True for
        members to leave out; a skipped folder (name ending ``/``) should also skip its contents.
    :param rename: function of a member name returning the name to extract it as.
    """

    def __init__(self, extract_dir, skip=None, rename=None):
        super().__init__(extract_dir=extract_dir)
        self.skip = skip or (lambda name: False)
        self.rename = rename or (lambda name: name)

    def _kept(self, fname):
        with _zipfile.ZipFile(fname, 'r') as archive:
            return [m for m in archive.namelist() if not self.skip(m)]

    def _all_members(self, fname):
        # what pooch checks for in the extract folder before deciding to extract again
        return [self.rename(m) for m in self._kept(fname)]

    def _extract_file(self, fname, extract_dir):
        with _zipfile.ZipFile(fname, 'r') as archive:
            for member in archive.namelist():
                if self.skip(member):
                    continue
                target = _os.path.join(extract_dir, self.rename(member))
                if member.endswith('/'):
                    _os.makedirs(target, exist_ok=True)
                    continue
                _os.makedirs(_os.path.dirname(target), exist_ok=True)
                with archive.open(member) as src, open(target, 'wb') as dst:
                    _shutil.copyfileobj(src, dst)


def skip_folder(name):
    """A ``skip`` for :class:`UnzipTidy`: leave out every member under a folder called ``name``."""
    return lambda member: ('/' + name + '/') in ('/' + member)


def underscores_in_file_names(member):
    """A ``rename`` for :class:`UnzipTidy`: spaces in the file name (not its folders) become
    underscores, which helps pygmt."""
    folder, _, base = member.rpartition('/')
    return (folder + '/' if folder else '') + base.replace(' ', '_')
