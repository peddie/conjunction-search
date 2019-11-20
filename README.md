conjunction-search
==================

Searches for astronomical conjunctions.

Requirements
------------

You will have to get a bunch of SPICE kernels from JPL:

```
$ grep bsp search.py
        self.ephemeris = load('de430.bsp')
    planet_barycenters = load('de430.bsp')
        SolarSystemBody('Mars', load('mar097.bsp')['mars'], ephem.Mars()),
        SolarSystemBody('Jupiter', load('jup310.bsp')['jupiter'], ephem.Jupiter()),
        SolarSystemBody('Saturn', load('sat319.bsp')['saturn'], ephem.Saturn()),
        SolarSystemBody('Uranus', load('ura111.bsp')['uranus'], ephem.Uranus()),
        SolarSystemBody('Neptune', load('nep081.bsp')['neptune'], ephem.Neptune()),
        SolarSystemBody('Pluto', load('plu055.bsp')['pluto'], ephem.Pluto()),

```

Then go ahead and do a search (takes a WHILE):

    python3 search.py

Find out how to control the search:

    python3 search.py --help