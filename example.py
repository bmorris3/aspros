import astropy.units as u
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt

from aspros import simulate_lc, inject_transits, bls_peakfinder

period = (3 * np.random.rand() + 3) * u.hour
epoch = Time('2020-04-01') + np.random.rand()*u.day
radius = 1500 * u.km
inc = 90 * u.deg
periods = np.linspace(2, 8, 1000) * u.hour
transit_duration = 2 * u.min

seed = None  # 42
clean_lc = simulate_lc(24*u.hour, efficiency=0.6, seed=seed)

transit_lc = inject_transits(clean_lc, period, epoch, radius, inc)

results, bests, stats = transit_lc.bls(periods=periods,
                                       duration=transit_duration)
best_period, best_duration, best_epoch = bests

# Plot results
fig, ax = plt.subplots(1, 3, figsize=(10, 4))
transit_lc.plot(ax=ax[0])

top_powers, significance = bls_peakfinder(results)

ax[1].plot(results.period.to(u.hour), results.power)
ax[1].axvline(results.period.to(u.hour)[top_powers[0]].value,
              ls=':', color='k', zorder=-10)
ax[1].axvline(results.period.to(u.hour)[top_powers[1]].value,
              ls=':', color='r', zorder=-12)

ax[1].set_xlabel('Period [hour]')
ax[1].set_ylabel('BLS Power')

phases = (((transit_lc.times.jd - best_epoch.jd) %
           best_period.to(u.day).value) / best_period.to(u.day).value)
phases[phases > 0.5] -= 1
ax[2].errorbar(phases, transit_lc.fluxes, transit_lc.errors, color='k',
               fmt='.', ecolor='silver')
ax[2].set_xlim([-0.02, 0.02])
ax[2].set(xlabel='Phase', ylabel='Flux')
ax[2].axvspan(-0.001, 0.001, color='k', alpha=0.1)

ax[0].set_title('Period = {0:.2f}; Radius = {1:.0f}'.format(period, radius))
ax[1].set_title('Peak height ratio = {0:.1f}'.format(significance))
# ax[2].set_title('Depth = {0:.2f}'
#                 .format(np.median(transit_lc.fluxes) -
#                         transit_lc.fluxes.min()))

intransit = np.abs(phases) < 0.0015
flux_intransit = np.median(transit_lc.fluxes[intransit])
flux_oot = np.median(transit_lc.fluxes)

snr = ((flux_oot - flux_intransit) / (np.median(transit_lc.errors) /
                                      np.count_nonzero(intransit)**0.5))

ax[2].set_title('S/N = {0:.1f}'.format(snr))

fig.tight_layout()
plt.show()
