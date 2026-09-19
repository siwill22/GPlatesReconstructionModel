"""Liveness of the dataset URLs.

Nothing here tests gprm's own logic. It exists because the dataset loaders point at 50-odd
fixed URLs, several of them on personal academic domains, and the failure mode when one rots
is a raw traceback from deep inside pooch. Running this on a schedule turns that into an
early warning.

Marked ``network`` and excluded from the default run.
"""
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.network

DATASETS_DIR = Path(__file__).parent.parent / 'gprm' / 'datasets'
URL_PATTERN = re.compile(r'url\s*=\s*["\'](https?://[^"\']+)["\']')


def collect_urls():
    """Every url= argument in the dataset modules, with the file and line it came from."""
    found = []
    for path in sorted(DATASETS_DIR.glob('*.py')):
        for lineno, line in enumerate(path.read_text().splitlines(), start=1):
            if line.lstrip().startswith('#'):
                continue        # commented-out URLs are not live dependencies
            match = URL_PATTERN.search(line)
            if match:
                found.append(pytest.param(match.group(1),
                                          id='{}:{}'.format(path.name, lineno)))
    return found


@pytest.mark.parametrize('url', collect_urls())
def test_dataset_url_is_reachable(url):
    requests = pytest.importorskip('requests')

    try:
        response = requests.head(url, allow_redirects=True, timeout=30)
        # Some servers refuse HEAD but serve GET perfectly well
        if response.status_code >= 400:
            response = requests.get(url, stream=True, allow_redirects=True, timeout=30)
            response.close()
    except requests.RequestException as error:
        pytest.fail('{} is unreachable: {}'.format(url, error))

    if response.status_code in (401, 403, 429):
        # The host answered but declined an automated request. Some repositories (Dryad, for
        # one) sit behind bot protection and still serve pooch perfectly well, so this is
        # inconclusive rather than a rotted link.
        pytest.skip('{} returned HTTP {} to an automated request'.format(url,
                                                                        response.status_code))

    assert response.status_code < 400, '{} returned HTTP {}'.format(url, response.status_code)


def test_no_new_plaintext_http_urls():
    """Plain http:// URLs are the ones most likely to rot, and are worth not adding more of.

    The existing ones are grandfathered in by this count; if it goes up, something new was
    added over http when it could have used https.
    """
    plaintext = [param.values[0] for param in collect_urls()
                 if param.values[0].startswith('http://')]

    assert len(plaintext) <= 7, 'new plaintext http:// dataset URLs added: {}'.format(plaintext)
