from batman import TransitModel, TransitParams
import numpy as np
import astropy.units as u
from astropy.constants import R_earth, G, M_sun

from .lightcurve import LightCurve

__all__ = ['inject_transits', 'period_to_a']

M_wd = 0.6 * M_sun
R_wd = 1 * R_earth


@u.quantity_input(radius=u.m, period=u.d, inc=u.deg, rstar=u.m)
def inject_transits(lc, period, epoch, radius, inc, rstar=1*R_earth, a=None):
    """

    Parameters
    ----------
    lc : `~aspros.LightCurve`
    period : `~astropy.units.Quantity`
    epoch : `~astropy.time.Time`
    radius : `~astropy.units.Quantity`
    inc : `~astropy.units.Quantity`
    a : float (optional)
    rstar : `~astropy.units.Quantity` (optional)

    Returns
    -------
    lc_transit : `~aspros.LightCurve`
    """
    p = TransitParams()
    p.per = period.to(u.day).value
    p.rp = float(radius / rstar)
    p.ecc = 0
    p.w = 90
    p.inc = inc.to(u.deg).value

    if a is None:
        a = period_to_a(period)

    p.a = a
    p.t0 = epoch.jd
    p.u = [0.4, 0.2]
    p.limb_dark = 'quadratic'

    model = TransitModel(p, lc.times.jd, supersample_factor=10,
                         exp_time=(60*u.s).to(u.day).value).light_curve(p)

    new_lc = LightCurve(lc.times, lc.fluxes * model, lc.errors)
    return new_lc


@u.quantity_input(period=u.d)
def period_to_a(period):
    a_rs = (period**2 * (G * M_wd) / 4 / np.pi**2)**(1/3) / R_wd
    return a_rs
