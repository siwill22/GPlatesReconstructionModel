"""UnzipTidy: archives are tidied while they are extracted, so a later load extracts nothing.

Before it, fetch_Muller2019 and fetch_Clennett deleted DeformingMeshPoints after extracting and
fetch_Paleomap deleted .cache files and renamed files with spaces. pooch extracts an archive
again whenever a member is missing, so every call unpacked the whole archive again and gave its
files new modification times.
"""
import os
import zipfile

from gprm.datasets._fetch import UnzipTidy, skip_folder, underscores_in_file_names

MEMBERS = ['model/', 'model/rotations.rot', 'model/DeformingMeshPoints/', 'model/DeformingMeshPoints/p.gpml',
           'model/grids/', 'model/grids/Map 1 PaleoDEM_10Ma.nc', 'model/grids/old.cache']


def _archive(tmp_path):
    path = tmp_path / 'model.zip'
    with zipfile.ZipFile(path, 'w') as z:
        for m in MEMBERS:
            if m.endswith('/'):
                z.writestr(m, '')
            else:
                z.writestr(m, 'contents of ' + m)
    return str(path)


def _processor():
    def skip(m):
        return skip_folder('DeformingMeshPoints')(m) or m.endswith('.cache')
    return UnzipTidy(extract_dir='extracted', skip=skip, rename=underscores_in_file_names)


def test_skipped_members_are_left_out_and_renamed_ones_arrive_renamed(tmp_path):
    fnames = _processor()(_archive(tmp_path), 'download', None)
    rel = sorted(os.path.relpath(f, tmp_path / 'extracted') for f in fnames)
    assert rel == ['model/grids/Map_1_PaleoDEM_10Ma.nc', 'model/rotations.rot']
    assert not (tmp_path / 'extracted' / 'model' / 'DeformingMeshPoints').exists()


def test_a_second_load_extracts_nothing(tmp_path):
    archive = _archive(tmp_path)
    _processor()(archive, 'download', None)
    target = tmp_path / 'extracted' / 'model' / 'rotations.rot'
    os.utime(target, (1_000_000_000, 1_000_000_000))       # mark it; a re-extraction would reset it
    calls = []
    p = _processor()
    p._extract_file = lambda *a: calls.append(a)
    p(archive, 'fetch', None)
    assert calls == []
    assert os.stat(target).st_mtime == 1_000_000_000


def test_a_missing_kept_member_is_extracted_again(tmp_path):
    archive = _archive(tmp_path)
    _processor()(archive, 'download', None)
    os.remove(tmp_path / 'extracted' / 'model' / 'rotations.rot')
    _processor()(archive, 'fetch', None)
    assert (tmp_path / 'extracted' / 'model' / 'rotations.rot').exists()


def test_rename_touches_the_file_name_only():
    assert underscores_in_file_names('a b/c d.nc') == 'a b/c_d.nc'
    assert underscores_in_file_names('c d.nc') == 'c_d.nc'
    assert skip_folder('X')('top/X/file') and skip_folder('X')('X/') and not skip_folder('X')('top/XY/file')
