import astropy.units as u
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt

from aspros import simulate_lc, inject_transits, bls_peakfinder

seed = None  # 42
clean_lc = simulate_lc(24*u.hour, efficiency=0.6, seed=seed)

period = (3 * np.random.rand() + 3) * u.hour
epoch = Time('2020-04-01') + np.random.rand()*u.day
radius = 1500 * u.km
inc = 90 * u.deg

transit_lc = inject_transits(clean_lc, period, epoch, radius, inc)

periods = np.linspace(2, 8, 1000) * u.hour
results = transit_lc.bls(periods=periods, duration=2*u.min)

fig, ax = plt.subplots(1, 2, figsize=(6, 3))
transit_lc.plot(ax=ax[0])

top_powers, significance = bls_peakfinder(results)

print('Period:', period)
print('Peak height ratio:', significance)

ax[1].plot(results.period.to(u.hour), results.power)

ax[1].axvline(results.period.to(u.hour)[top_powers[0]].value,
              ls=':', color='k', zorder=-10)
ax[1].axvline(results.period.to(u.hour)[top_powers[1]].value,
              ls=':', color='r', zorder=-12)

ax[1].set_xlabel('Period [hour]')
ax[1].set_ylabel('BLS Power')
fig.tight_layout()
plt.show()