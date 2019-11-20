conjunction-search
==================

Searches for astronomical conjunctions.

Requirements
------------

It's written in python 3.  You can get the needed libraries by

    pip3 install -r requirements.txt

You might have to either prefix this with `sudo` or also use the `--user` flag;
I haven't yet understood the right approach, if one exists.

You will have to get a bunch of SPICE kernels from JPL.  In theory the program
should fetch all this stuff itself when run, but this functionality has been
pretty unreliable for me -- consider yourself warned.

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

It will also fetch the Hipparcos catalog and some leap second and EOP data.

Running
-------

Go ahead and do a search (takes a WHILE):

    python3 search.py

Find out more, including how to control the search:

    python3 search.py --help
