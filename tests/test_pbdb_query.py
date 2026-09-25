"""pbdb_query reloads a saved snapshot without touching the network."""

from gprm.datasets.Strat import pbdb_query


def test_existing_snapshot_is_read_not_downloaded(tmp_path, monkeypatch):
    import requests

    def no_network(*args, **kwargs):
        raise AssertionError('pbdb_query downloaded although the snapshot exists')

    monkeypatch.setattr(requests, 'get', no_network)
    snapshot = tmp_path / 'reefs.csv'
    snapshot.write_text('collection_no,lng,lat,environment,max_ma,min_ma\n'
                        '1,-3.5,50.3,"reef, buildup or bioherm",393.3,382.7\n'
                        '2,150.0,-20.0,perireef or subreef,5.3,2.6\n')

    gdf = pbdb_query('colls/list', {'envtype': 'reef'}, snapshot=str(snapshot))

    assert list(gdf['collection_no']) == [1, 2]
    assert list(gdf['Longitude']) == [-3.5, 150.0]
    assert gdf.attrs['pbdb_snapshot'] == str(snapshot)
    assert gdf.attrs['gprm_age']['range'] == ('min_ma', 'max_ma')
