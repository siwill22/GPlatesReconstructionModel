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
from pathlib import Path as _Path

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
