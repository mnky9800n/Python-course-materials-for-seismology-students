from cartopy.io.img_tiles import GoogleTiles


class ShadedReliefESRI(GoogleTiles):
    # shaded relief
    def _image_url(self, tile):
        x, y, z = tile
        url = ('https://server.arcgisonline.com/ArcGIS/rest/services/' \
               'World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg').format(
               z=z, y=y, x=x)
        return url

class notSure(GoogleTiles):
    #
    def _image_url(self, tile):
        x, y, z = tile
        url = ('https://gibs.earthdata.nasa.gov/wmts/epsg3857/best/MODIS_Terra_Aerosol/' \
               'default/2014-04-09/GoogleMapsCompatible_Level6/{z}/{y}/{x}.png').format(z=z, y=y, x=x)
        return url
