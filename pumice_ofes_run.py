from parcels import FieldSet, Field, ParticleSet, JITParticle, ErrorCode, AdvectionRK4, BrownianMotion2D
from datetime import timedelta as delta
from glob import glob
import numpy as np
import xlrd


def set_fields():
    ofesufiles = sorted(glob('/Volumes/data01/OFESdata/OFES_0.1_HIND/allveldata/nest_1_200[7-9]*u.nc'))
    dimensions = {'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time'}

    uofes = Field.from_netcdf(ofesufiles, 'zu', dimensions, fieldtype='U')

    ofesvfiles = [f.replace('u.nc', 'v.nc') for f in ofesufiles]
    vofes = Field.from_netcdf(ofesvfiles, 'zv', dimensions, fieldtype='V')

    stokesfiles = sorted(glob('/Volumes/data01/WaveWatch3data/WW3-GLOB-30M_200[7-9]*'))
    dimensions = {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'}
    uuss = Field.from_netcdf(stokesfiles, 'uuss', dimensions, fieldtype='U')
    vuss = Field.from_netcdf(stokesfiles, 'vuss', dimensions, fieldtype='V')

    # landpfiles = sorted(glob('/Volumes/data01/OFESdata/OFES_auxiliaryfiles/boundary_velocities.nc'))
    # uland = Field.from_netcdf(landpfiles, 'MaskUvel', dimensions, allow_time_extrapolation=True, fieldtype='U')
    # vland = Field.from_netcdf(landpfiles, 'MaskVvel', dimensions, allow_time_extrapolation=True, fieldtype='V')
    # fieldset = FieldSet(U=[uofes, uuss, uland], V=[vofes, vuss, vland])
    fieldset = FieldSet(U=[uofes, uuss], V=[vofes, vuss])

    fieldset.add_periodic_halo(zonal=True, meridional=False, halosize=5)
    return fieldset


def WrapLon(particle, fieldset, time, dt):
    if particle.lon > 360.:
        particle.lon = particle.lon - 360.
    if particle.lon < 0.:
        particle.lon = particle.lon + 360.


def OutOfBounds(particle, fieldset, time, dt):
    particle.delete()


fieldset = set_fields()

size2D = (fieldset.U[0].grid.ydim, fieldset.U[0].grid.xdim)
fieldset.add_field(Field('Kh_zonal', data=10 * np.ones(size2D), lon=fieldset.U[0].grid.lon, lat=fieldset.U[0].grid.lat,
                         mesh='spherical', allow_time_extrapolation=True))
fieldset.add_field(Field('Kh_meridional', data=10 * np.ones(size2D), lon=fieldset.U[0].grid.lon, lat=fieldset.U[0].grid.lat,
                         mesh='spherical', allow_time_extrapolation=True))

nperloc = 100
rundays = 729

book = xlrd.open_workbook("pumicesources.xlsx")
sh = book.sheet_by_index(0)
lons = [sh.cell_value(rowx=rx, colx=3) for rx in range(sh.nrows)]
lats = [sh.cell_value(rowx=rx, colx=2) for rx in range(sh.nrows)]
lons = [ln if ln > 0 else ln + 360 for ln in lons]

pset = ParticleSet.from_list(fieldset=fieldset, pclass=JITParticle, lon=np.tile(lons, [nperloc]),
                             lat=np.tile(lats, [nperloc]))

ofile = pset.ParticleFile(name='pumice_wstokes', outputdt=delta(days=5))

kernels = pset.Kernel(WrapLon) + AdvectionRK4 + pset.Kernel(BrownianMotion2D)
pset.execute(kernels, runtime=delta(days=rundays), dt=delta(hours=1),
             output_file=ofile, recovery={ErrorCode.ErrorOutOfBounds: OutOfBounds})
