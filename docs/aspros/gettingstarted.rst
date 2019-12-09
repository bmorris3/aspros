Getting Started
===============

Let's start by importing the necessary packages

.. code-block:: python

    from aspros import simulate_lc, inject_transits

    import astropy.units as u
    from astropy.time import Time

    import numpy as np
    import matplotlib.pyplot as plt

Simulate a light curve with CHEOPS-like noise and window function:

.. code-block:: python

    seed = 42  # Makes random number generator reproducible
    duration = 24 * u.hour
    efficiency = 0.6
    clean_lc = simulate_lc(duration, efficiency=efficiency, seed=seed)

Inject a transiting object with specified orbital properties:

.. code-block:: python

    period = 3 * u.hour
    epoch = Time('2020-04-01') + np.random.rand()*u.day
    radius = 1500 * u.km
    inc = 90 * u.deg

    transit_lc = inject_transits(clean_lc, period, epoch, radius, inc)

Construct a Box Least Squares periodogram and inspect it for peaks:

.. code-block:: python

    periods = np.linspace(2, 12, 1500) * u.hour
    results, bests, stats = transit_lc.bls(periods=periods, duration=2*u.min)

    fig, ax = plt.subplots(1, 2, figsize=(8, 2))
    transit_lc.plot(ax=ax[0])
    ax[1].plot(results.period.to(u.hour), results.power)
    ax[1].set_xlabel('Period [hour]')
    ax[1].set_ylabel('BLS Power')
    fig.tight_layout()
    plt.show()

.. plot::

    import astropy.units as u
    from astropy.time import Time
    from aspros import simulate_lc, inject_transits
    import numpy as np
    import matplotlib.pyplot as plt

    seed = 42
    clean_lc = simulate_lc(24*u.hour, efficiency=0.6, seed=seed)

    period = 3 * u.hour
    epoch = Time('2020-04-01') + np.random.rand()*u.day
    radius = 1500 * u.km
    inc = 90 * u.deg

    transit_lc = inject_transits(clean_lc, period, epoch, radius, inc)
    periods = np.linspace(2, 12, 1500) * u.hour
    results, bests, stats = transit_lc.bls(periods=periods, duration=2*u.min)

    fig, ax = plt.subplots(1, 2, figsize=(8, 2))
    transit_lc.plot(ax=ax[0])
    ax[1].plot(results.period.to(u.hour), results.power)
    ax[1].set_xlabel('Period [hour]')
    ax[1].set_ylabel('BLS Power')
    fig.tight_layout()
    plt.show()
