import numpy as np
import os

import astropy.units as u
from astropy.time import Time

from .lightcurve import LightCurve

__all__ = ['simulate_lc']


start_time = Time('2020-04-01')
cheops_period = 100 * u.min


@u.quantity_input(duration=u.s)
def simulate_lc(duration, efficiency=0.5, seed=None, noise_scale_factor=2):
    """
    Simulate a white dwarf light curve observed with CHEOPS.

    Parameters
    ----------
    duration : `~astropy.units.Quantity`
        Duration of the simulated observations
    efficiency : float (optional)
        Efficiency of the observations, defaults to 0.5.
    seed : int (optional)
        Random seed
    noise_scale_factor : float (optional)
        Scale up the observed noise by a factor of the estimated photon noise,
        defaults to 2.

    Returns
    -------
    lc : `~aspros.LightCurve`
        Light curve of the object
    """

    if seed is not None:
        np.random.seed(seed)

    path = os.path.join(os.path.dirname(__file__), 'data', 'count_rate.txt')
    count_rate = float(open(path, 'r').read()) * u.ct / u.s
    signal = (count_rate * 60*u.s).to(u.ct).value
    noise = noise_scale_factor / np.sqrt(signal)

    times = start_time + np.arange(0, duration.to(u.s).value, 60) * u.s
    fluxes = noise * np.random.randn(len(times)) + 1
    errors = noise * np.ones(len(times))

    # Earth occultation mask
    epoch = start_time + np.random.rand() * u.day
    phases = (((times.jd - epoch.jd) % cheops_period.to(u.day).value) /
              cheops_period.to(u.day).value)

    earth_mask = phases < efficiency

    return LightCurve(times[earth_mask], fluxes[earth_mask],
                      errors[earth_mask])
