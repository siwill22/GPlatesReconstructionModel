"""Fetch one member of a remote ZIP archive without downloading the whole archive.

Some published supplements bundle everything -- data, model output, notebooks, figures --
into a single very large zip. The Flament et al. (2022) supplement on Zenodo is 772 MB, and
the two eruption catalogues gprm wants out of it total about 7 MB. Running that through
pooch's ``Unzip`` processor, as the other loaders do, would cost every user the full download
plus well over a gigabyte of unpacked cache.

A zip's central directory lives at the end of the file, and each member is stored
independently, so with HTTP range requests it is possible to read the directory and then
inflate just the member wanted. That is all this module does: a seekable file object backed
by range requests (stdlib only -- no new dependency), handed to :mod:`zipfile`.

This is only worth it for a big archive holding a small file. For everything else use
``_fetch.retrieve`` with pooch's ``Unzip``, which stays the default idiom in this package.

MIT License

Copyright (c) 2017-2021 Simon Williams
"""

import hashlib as _hashlib
import io as _io
import os as _os
import urllib.request as _urllib_request
import zipfile as _zipfile
from pathlib import Path as _Path

from ._fetch import fetch_error as _fetch_error

_USER_AGENT = 'gprm (https://github.com/siwill22/GPlatesReconstructionModel)'
_TIMEOUT = 60


class RangeRequestsUnsupported(RuntimeError):
    """The server will not serve byte ranges, so a member cannot be extracted on its own."""


def _urlopen(url, headers=None):
    request = _urllib_request.Request(url, headers=dict({'User-Agent': _USER_AGENT},
                                                        **(headers or {})))
    return _urllib_request.urlopen(request, timeout=_TIMEOUT)


def _content_length(url):
    request = _urllib_request.Request(url, headers={'User-Agent': _USER_AGENT}, method='HEAD')
    with _urllib_request.urlopen(request, timeout=_TIMEOUT) as response:
        length = response.headers.get('Content-Length')
    if length is None:
        raise RangeRequestsUnsupported(
            'The server did not report a Content-Length for {}, so the end of the zip '
            '(where its index lives) cannot be located.'.format(url))
    return int(length)


def _require_range_support(url):
    """Ask for a single byte. A server that ignores Range answers 200 with the whole file.

    Checked up front, and the body deliberately not read, so that an unsupported server costs
    one aborted connection rather than a silent multi-hundred-megabyte download.
    """
    with _urlopen(url, {'Range': 'bytes=0-0'}) as response:
        status = getattr(response, 'status', None) or response.getcode()
    if status != 206:
        raise RangeRequestsUnsupported(
            'The server answered a byte-range request with HTTP {} instead of 206, so it is '
            'serving the whole file rather than the requested range. gprm will not download '
            'the entire archive to reach one file inside it.'.format(status))


class _HttpRangeFile(_io.RawIOBase):
    """A read-only, seekable view of a remote file, backed by HTTP range requests."""

    def __init__(self, url, size=None):
        self.url = url
        self.size = _content_length(url) if size is None else size
        self._pos = 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self._pos

    def seek(self, offset, whence=_os.SEEK_SET):
        if whence == _os.SEEK_SET:
            new = offset
        elif whence == _os.SEEK_CUR:
            new = self._pos + offset
        elif whence == _os.SEEK_END:
            new = self.size + offset
        else:
            raise ValueError('Unsupported whence value {!r}'.format(whence))
        self._pos = max(0, min(new, self.size))
        return self._pos

    def readinto(self, buffer):
        wanted = len(buffer)
        if wanted == 0 or self._pos >= self.size:
            return 0
        last = min(self._pos + wanted, self.size) - 1
        with _urlopen(self.url, {'Range': 'bytes={}-{}'.format(self._pos, last)}) as response:
            data = response.read()
        buffer[:len(data)] = data
        self._pos += len(data)
        return len(data)


def _split_hash(known_hash):
    """Accept pooch's 'alg:hex' spelling, and a bare hex digest meaning sha256."""
    if ':' in known_hash:
        algorithm, _, digest = known_hash.partition(':')
        return algorithm.lower(), digest.lower()
    return 'sha256', known_hash.lower()


def _digest(data, algorithm):
    hasher = _hashlib.new(algorithm)
    hasher.update(data)
    return hasher.hexdigest()


def _cached_digest(target, algorithm, chunk=1 << 20):
    hasher = _hashlib.new(algorithm)
    with open(target, 'rb') as handle:
        for block in iter(lambda: handle.read(chunk), b''):
            hasher.update(block)
    return hasher.hexdigest()


def retrieve_zip_member(url, member, known_hash, path, fname=None):
    """Download one member of a remote zip and cache it, as ``_fetch.retrieve`` caches a file.

    :param url: the zip archive's URL. The server must honour byte ranges.
    :param member: the member's full path inside the archive.
    :param known_hash: checksum of the **extracted member**, as ``'sha256:...'``, ``'md5:...'``
        or a bare sha256 digest. Not the checksum of the archive.
    :param path: cache directory, normally ``pooch.os_cache('gprm')``.
    :param fname: name to cache the member under. Defaults to its basename inside the archive,
        which must therefore be distinctive enough not to collide with another dataset.
    :returns: str path to the cached file.
    :raises gprm.datasets.DatasetFetchError: if the archive cannot be read, the member is not
        in it, or the extracted bytes do not match ``known_hash``.
    """
    algorithm, expected = _split_hash(known_hash)
    directory = _Path(path)
    target = directory / (member.rsplit('/', 1)[-1] if fname is None else fname)

    if target.exists() and _cached_digest(target, algorithm) == expected:
        return str(target)

    try:
        _require_range_support(url)
        archive = _zipfile.ZipFile(_io.BufferedReader(_HttpRangeFile(url), buffer_size=1 << 22))
        data = archive.read(member)
    except Exception as err:
        raise _fetch_error(err, url, path=directory,
                           extra=['member: {}'.format(member)]) from err

    actual = _digest(data, algorithm)
    if actual != expected:
        raise _fetch_error(
            ValueError('{} hash of the extracted member is {}, expected {}'.format(
                algorithm, actual, expected)),
            url, path=directory, extra=['member: {}'.format(member)])

    directory.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + '.part')
    partial.write_bytes(data)
    _os.replace(partial, target)      # atomic, so an interrupted run leaves no half file
    return str(target)
