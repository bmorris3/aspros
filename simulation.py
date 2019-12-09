import astropy.units as u
from astropy.time import Time
import numpy as np

from aspros import simulate_lc, inject_transits, bls_peakfinder

n_trials = 1000000

detected = []
not_detected = []

for i in range(n_trials):
    period = (9 * np.random.rand() + 3) * u.hour
    epoch = Time('2020-04-01') + np.random.rand() * u.day
    radius = (500 + 2500 * np.random.rand()) * u.km
    inc = 90 * u.deg
    periods = np.linspace(2, 12, 1500) * u.hour
    transit_duration = 2 * u.min

    seed = None  # 42
    clean_lc = simulate_lc(24*u.hour, efficiency=0.6, seed=seed)

    transit_lc = inject_transits(clean_lc, period, epoch, radius, inc)

    results, bests, stats = transit_lc.bls(periods=periods,
                                           duration=transit_duration)
    best_period, best_duration, best_epoch = bests

    top_powers, significance = bls_peakfinder(results)

    # Compute S/N
    phases = (((transit_lc.times.jd - best_epoch.jd) %
               best_period.to(u.day).value) / best_period.to(u.day).value)
    phases[phases > 0.5] -= 1
    intransit = np.abs(phases) < 0.0015
    flux_intransit = np.median(transit_lc.fluxes[intransit])
    flux_oot = np.median(transit_lc.fluxes)
    snr = ((flux_oot - flux_intransit) / (np.median(transit_lc.errors) /
                                          np.count_nonzero(intransit) ** 0.5))

    # if significance > 1:
    if snr > 5:
        detected.append([period.to(u.d).value, radius.to(u.km).value, snr])
    else:
        not_detected.append([period.to(u.d).value, radius.to(u.km).value, snr])

np.save('outputs/detected.npy', detected)
np.save('outputs/not_detected.npy', not_detected)
