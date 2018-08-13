from parcels import (FieldSet, Field, ParticleSet, JITParticle, ErrorCode,
                     AdvectionRK4, BrownianMotion2D, Variable)
from datetime import timedelta as delta
from glob import glob
import math
import numpy as np
import xlrd


def set_fields():
    datestr = '201[0-2]'
    hycomfiles = sorted(glob('/Volumes/data01/HYCOMdata/GLBu0.08_expt_19.1_surf/hycom_GLBu0.08_191_%s*' % datestr))
    dimensions = {'lat': 'lat', 'lon': 'lon', 'time': 'time'}
    uhycom = Field.from_netcdf(hycomfiles, 'water_u', dimensions, fieldtype='U')
    vhycom = Field.from_netcdf(hycomfiles, 'water_v', dimensions, fieldtype='V')
    uhycom.vmin = -99.
    uhycom.set_scaling_factor(0.001)
    vhycom.set_scaling_factor(0.001)
    vhycom.vmin = -99.

    stokesfiles = sorted(glob('/Volumes/data01/WaveWatch3data/WW3-GLOB-30M_%s*' % datestr))
    dimensions = {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'}
    uuss = Field.from_netcdf(stokesfiles, 'uuss', dimensions, fieldtype='U')
    vuss = Field.from_netcdf(stokesfiles, 'vuss', dimensions, fieldtype='V')

    fieldset = FieldSet(U=[uhycom, uuss], V=[vhycom, vuss])

    fieldset.add_periodic_halo(zonal=True, meridional=False, halosize=5)
    return fieldset


def WrapLon(particle, fieldset, time, dt):
    if particle.lon > 360.:
        particle.lon = particle.lon - 360.
    if particle.lon < 0.:
        particle.lon = particle.lon + 360.


def OutOfBounds(particle, fieldset, time, dt):
    particle.delete()


def Age(particle, fieldset, time, dt):
    particle.age = particle.age + math.fabs(dt)
    if particle.age > fieldset.maxage:
        particle.delete()


fieldset = set_fields()
fieldset.add_constant('maxage', 730.*86400)

size2D = (fieldset.U[0].grid.ydim, fieldset.U[0].grid.xdim)
fieldset.add_field(Field('Kh_zonal', data=10 * np.ones(size2D), lon=fieldset.U[0].grid.lon, lat=fieldset.U[0].grid.lat,
                         mesh='spherical', allow_time_extrapolation=True))
fieldset.add_field(Field('Kh_meridional', data=10 * np.ones(size2D), lon=fieldset.U[0].grid.lon, lat=fieldset.U[0].grid.lat,
                         mesh='spherical', allow_time_extrapolation=True))

book = xlrd.open_workbook("pumicesources.xlsx")
sh = book.sheet_by_index(0)
lons = [sh.cell_value(rowx=rx, colx=3) for rx in range(sh.nrows)]
lats = [sh.cell_value(rowx=rx, colx=2) for rx in range(sh.nrows)]
lons = [ln if ln > 0 else ln + 360 for ln in lons]
times = fieldset.U[0].time[0]+np.arange(0, 365*86400, 86400*5)


class PumiceParticle(JITParticle):
    age = Variable('age', dtype=np.float32, initial=0.)


pset = ParticleSet.from_list(fieldset=fieldset, pclass=PumiceParticle, lon=np.tile(lons, [len(times)]),
                             lat=np.tile(lats, [len(times)]), time=np.repeat(times, len(lons)))

ofile = pset.ParticleFile(name='pumice_wstokes_hycom_delayedtime', outputdt=delta(days=5))

kernels = pset.Kernel(WrapLon) + AdvectionRK4 + BrownianMotion2D + Age
pset.execute(kernels, endtime=fieldset.U[0].time[-1], dt=delta(hours=1),
             output_file=ofile, recovery={ErrorCode.ErrorOutOfBounds: OutOfBounds})
