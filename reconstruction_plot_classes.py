import os
import pygplates
from ptt.resolve_topologies import resolve_topologies as topology2gmt
from ptt.utils.call_system_command import call_system_command


class gmt_reconstruction(object):

    def __init__(self, reconstruction_model, output_dir, output_file_stem):
        self.reconstruction_model = reconstruction_model
        self.output_dir = output_dir
        self.output_file_stem = output_file_stem

    def set_projection(self, projection):
        self.projection = projection

    def set_region(self, region):
        self.region = region

    def plot_snapshot(self, reconstruction_time, anchor_plate_id = 0, 
                      layers=['continents','coastlines','dynamic_polygons']):

        region = '%f/%f/%f/%f' % (self.region[0],self.region[1],self.region[2],self.region[3])

        outfile='%s/%s_%dMa.ps' % (self.output_dir,self.output_file_stem,reconstruction_time)

        call_system_command(['gmt','gmtset','COLOR_MODEL','RGB',
                                            'MAP_FRAME_TYPE','inside',
                                            'MAP_FRAME_PEN','1.0',
                                            'MAP_TICK_PEN_PRIMARY','1.0',
                                            'FONT_ANNOT_PRIMARY','9',
                                            'FONT_LABEL','8',
                                            'FONT_TITLE','8',
                                            'FONT_ANNOT_PRIMARY','Helvetica',
                                            'FORMAT_GEO_MAP','ddd'])

        call_system_command(['gmt', 'psbasemap', '-R%s' % region, self.projection,
                             '-Ba30f30::wesn', '-K', '>', outfile])
        
        if 'continents' in layers:
            output_reconstructed_continents_filename = 'tmp/continents.gmt'
            pygplates.reconstruct(self.reconstruction_model.continent_polygons, 
                                  self.reconstruction_model.rotation_model, 
                                  output_reconstructed_continents_filename, 
                                  reconstruction_time, anchor_plate_id=anchor_plate_id)
            call_system_command(['gmt', 'psxy', '-R%s' % region, self.projection, 
                                '-W0.1p,wheat', '-Gwheat', 'tmp/continents.gmt', '-O', '-K', '-N', '>>', outfile])

        if 'coastlines' in layers:
            output_reconstructed_coastlines_filename = 'tmp/coastlines.gmt'
            pygplates.reconstruct(self.reconstruction_model.coastlines, 
                                  self.reconstruction_model.rotation_model, 
                                  output_reconstructed_coastlines_filename, 
                                  reconstruction_time, anchor_plate_id=anchor_plate_id)
            call_system_command(['gmt', 'psxy', '-R', self.projection,
                                '-W0.2p,darkolivegreen', '-Gdarkolivegreen','-O', '-K', '-m', 'tmp/coastlines.gmt', '-V', '>>', outfile])
    
        if 'dynamic_polygons' in layers:
            output_filename_prefix = 'tmp/'
            output_filename_extension = 'gmt'
            topology2gmt(self.reconstruction_model.rotation_model,
                 self.reconstruction_model.dynamic_polygons,
                 reconstruction_time,
                 output_filename_prefix,
                 output_filename_extension,
                 anchor_plate_id)

            call_system_command(['gmt', 'psxy', '-R', self.projection,
                                 '-W0.6p,gray70', '-K', '-O', '-m', 'tmp/boundary_polygons_%0.2fMa.gmt' % reconstruction_time, '-V', '>>', outfile])
            call_system_command(['gmt', 'psxy', '-R', self.projection,
                                 '-W0.6p,red', '-K', '-O', '-m', 'tmp/ridge_transform_boundaries_%0.2fMa.gmt' % reconstruction_time, '-V', '>>', outfile])

            #plot subduction zones
            call_system_command(['gmt', 'psxy', '-R', self.projection,
                                 '-W0.6p,black', '-Gblack', '-Sf15p/4plt', '-K', '-O', '-m', 'tmp/subduction_boundaries_sL_%0.2fMa.gmt' % reconstruction_time, '-V', '>>', outfile])
            call_system_command(['gmt', 'psxy', '-R', self.projection,
                                 '-W0.6p,black', '-Gblack', '-Sf15p/4prt', '-K', '-O', '-m', 'tmp/subduction_boundaries_sR_%0.2fMa.gmt' % reconstruction_time, '-V', '>>', outfile])

        call_system_command(['gmt', 'psbasemap', '-R%s' % region, self.projection, 
                             '-Ba30f30::wesn', '-O', '-K', '>>', outfile])
        call_system_command(['gmt', 'psclip', '-C', '-O', '>>', outfile])

        #convert ps into raster, -E set the resolution
        call_system_command(['gmt', 'ps2raster', outfile, '-A0.2c', '-E300', '-Tg', '-P'])

        self.image_postscript = outfile
        self.image_file = '%s.png' % outfile[:-3]



