import requests
import pygplates

default_url = 'https://gws.gplates.org'
default_model = 'MULLER2016'


def FeatureCollection(model=default_model, layer='rotations', url=default_url, proxy=''):
    """
    function to access gpml data files directly from the gplates-web-service data store 
    and load them into pygplates feature collections.
    :param model: 
    :param layer: 'rotations', 'coastlines', 'static_polygons', 'plate_polygons'
    :param url:
    :param proxy:
    """

    if layer=='rotations':
        requested_file = 'request.rot'
    else:
        requested_file = 'request.gpmlz'

    with open(requested_file, 'wb') as handle:
        r = requests.get('%s/model/get_model_layer/?model=%s&layer=%s' % (url,model,layer),
                         proxies={'http':'%s' % proxy},
                         stream=True)

        if r.status_code!=200:
            error_string = 'Remote request returned with message %s' % r.text
            raise Exception(error_string)
        else:
            for block in r.iter_content(1024):
                handle.write(block)

        if layer=='rotations':
            result = pygplates.RotationModel('request.rot')
        else:
            result = pygplates.FeatureCollection('request.gpmlz')

    return result


# fetch reconstructed coastlines
def gws_coastlines(recon_time, model=default_model, url=default_url, proxies={'http':''}):

    r = requests.get('%s/reconstruct/coastlines/?time=%0.2f' % (url,recon_time),
                     proxies=proxies)
    cs = json.loads(r.text)
    
    return cs

# fetch reconstructed topological plate polygons
def gws_plate_polygons(recon_time, model=default_model, url=default_url, proxies={'http':''}):

    r = requests.get('%s/topology/plate_polygons/?time=%0.2f' % (url,recon_time),
                     proxies=proxies)
    pp = json.loads(r.text)

    return pp

# fetch reconstructed static polygons
def gws_static_polygons(recon_time, model=default_model, url=default_url, proxies={'http':''}):

     r = requests.get('%s/reconstruct/static_polygons/?time=%0.2f' % (url,recon_time),
                      proxies=proxies)
     sp = json.loads(r.text)

     return sp

# fetch reconstructed motion paths
def gws_motion_path(recon_time, seedpoints, fixplate, movplate, time_min, time_max, time_step,
                    model=default_model, url=default_url, proxies={'http':''}):
    
    r = requests.get('%s/reconstruct/motion_path/?model=%s&time=%0.2f&seedpoints=%s&timespec=%s&fixplate=%d&movplate=%s' % \
                     (url,model,recon_time,seedpoints,'%s,%s,%s' % (time_min,time_max,time_step),fixplate,movplate),
                     proxies=proxies)

    motion_path = json.loads(r.text)

    return motion_path

# fetch reconstructed plate velocities
def gws_velocities(recon_time,polygon_type='static',
                   model=default_model, velocity_type='MagAzim', url=default_url, proxies={'http':''}):

    if polygon_type=='topology':
        r = requests.get('%s/velocity/plate_polygons/?model=%s&time=%0.2f&velocity_type=%s' % \
                         (url,model,recon_time,velocity_type),
                         proxies=proxies)
    else:
        r = requests.get('%s/velocity/static_polygons/?model=%s&time=%0.2f&velocity_type=%s' % \
                         (url,model,recon_time,velocity_type),
                         proxies=proxies)

    velocities = json.loads(r.text)

    return velocities

def gws_model_list(url=default_url, proxies={'http':''}):

    r = requests.get('%s/list/models' % url, proxies=proxies)

    model_info = json.loads(r.text)

    return model_info


