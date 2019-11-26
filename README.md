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

Here's what it says at the time of writing.

```
usage: search.py [-h] [-a LAT] [-o LON] [--site-name NAME]
                 [--start-date START] [--num-nights NIGHTS | --end-date END]
                 [--max-star-magnitude MAXSTARMAG] [--max-magnitude MAXMAG]
                 [--max-magnitude-delta DELTA] [--max-angle MAXANGLE]
                 [--serial | --max-processes N]

Search for astronomical conjunctions.

By default:
 - location is Brisbane, QLD, Australia.
 - start date is tonight (2019-11-26).
 - end date is one week from now (2019-12-03).

Other default values are displayed below in parentheses.

Please note that the progress output is a bit wonky when manually setting
max_workers values greater than one and fewer than the number of cores.  We
apologize for the inconvenience.

optional arguments:
  -h, --help            show this help message and exit

OBSERVER CONFIGURATION:
  -a LAT, --latitude LAT
                        Latitude of the observation location ("27.49665 S")
  -o LON, --longitude LON
                        Longitude of the observation location ("152.9883 E")
  --site-name NAME      Name of the observing site (just makes output prettier) (Home)

SEARCH CONFIGURATION:
  --start-date START    Date of the first night (sunset) to search (2019-11-26)
  --num-nights NIGHTS   Number of nights after the start date to search (7)
  --end-date END        Date of the last sunrise up to which to search (calculated from --num-nights)
  --max-star-magnitude MAXSTARMAG
                        Maximum star magnitude to consider when searching (4.0)
  --max-magnitude MAXMAG
                        Maximum planet/moon magnitude to consider when searching (5.0)
  --max-magnitude-delta DELTA
                        Maximum magnitude discrepancy (between objects in a conjunction) to consider when searching (3.0)
  --max-angle MAXANGLE  Maximum planet/moon magnitude to consider when searching in degrees (1.0)
  --serial              Search nights serially (default is to do it in parallel)
  --max-processes N     Maximum number of parallel searchers to run (default is the number of cores)
```

Example usage
-------------

    ./search.py --start-date 2020-01-01 --end-date 2020-04-01 --site-name 'Brisbane'

You will then see some startup messages, followed by some progress bars which
get more disorganized the more cores your computer is using.  So far parallelism
in Python has been quite painful, and I'm still learning how to deal with it --
sorry for the confusing output.

Once the search finishes, you will see some beautiful ASCII like the following:

```
[   SEARCH] Mean scan time was 559.278s per night.
Displaying events in Brisbane, at 27.49665 S, 152.9883 E (Earth) between
    2020-01-01 00:00:00
 and
    2020-04-01 00:00:00
 :

2020-01-07 00:00:00+10:00 (2 events)
Conjunction between Mars and 78820:
                time: 2020-01-08 04:47:48.354000+10:00
            distance: 0.874 degrees
magnitude difference: 1.020 (2.185 times brighter)
   maximum magnitude: 2.560
               score: 3.59525

Conjunction between Mars and 78933:
                time: 2020-01-08 04:47:48.354000+10:00
            distance: 0.988 degrees
magnitude difference: 2.390 (6.244 times brighter)
   maximum magnitude: 3.930
               score: 6.33725

2020-01-08 00:00:00+10:00 (2 events)
Conjunction between Mars and 78820:
                time: 2020-01-09 02:17:54.646000+10:00
            distance: 0.727 degrees
magnitude difference: 1.030 (2.202 times brighter)
   maximum magnitude: 2.560
               score: 3.60270

Conjunction between Mars and 78933:
                time: 2020-01-09 04:47:54.646000+10:00
            distance: 0.317 degrees
magnitude difference: 2.400 (6.292 times brighter)
   maximum magnitude: 3.930
               score: 6.33554

2020-01-09 00:00:00+10:00 (1 event)
Conjunction between Mars and 78933:
                time: 2020-01-10 02:17:59.361000+10:00
            distance: 0.298 degrees
magnitude difference: 2.400 (6.292 times brighter)
   maximum magnitude: 3.930
               score: 6.33520

2020-01-10 00:00:00+10:00 (1 event)
Conjunction between Mars and 78933:
                time: 2020-01-11 02:18:02.495000+10:00
            distance: 0.969 degrees
magnitude difference: 2.410 (6.341 times brighter)
   maximum magnitude: 3.930
               score: 6.35692

2020-03-19 00:00:00+10:00 (1 event)
Conjunction between Mars and Jupiter:
                time: 2020-03-20 05:30:10.017000+10:00
            distance: 0.780 degrees
magnitude difference: 2.840 (8.816 times brighter)
   maximum magnitude: 0.920
               score: 3.77361

2020-03-20 00:00:00+10:00 (1 event)
Conjunction between Mars and Jupiter:
                time: 2020-03-21 00:59:03.016000+10:00
            distance: 0.708 degrees
magnitude difference: 2.830 (8.749 times brighter)
   maximum magnitude: 0.910
               score: 3.75236

2020-03-21 00:00:00+10:00 (1 event)
Conjunction between Mars and Jupiter:
                time: 2020-03-22 00:57:55.946000+10:00
            distance: 0.952 degrees
magnitude difference: 2.830 (8.749 times brighter)
   maximum magnitude: 0.900
               score: 3.74661

2020-03-31 00:00:00+10:00 (1 event)
Conjunction between Mars and Saturn:
                time: 2020-04-01 03:16:46.980000+10:00
            distance: 0.902 degrees
magnitude difference: 0.120 (1.096 times brighter)
   maximum magnitude: 0.780
               score: 0.91575
```

This is intended to be enough information to get you started.  I don't provide
ra / dec; my assumption is that you will probably want to use a dedicated
program like xephem or stellarium to plan your astrophotography session once you
know what bodies to look for when.

Interpreting the output
-----------------------

If you find a conjunction with a star, your output will contain a Hipparcos ID.
ESA has a catalog search [available here](https://www.cosmos.esa.int/web/hipparcos/search-facility);
currently I don't do anything to make these IDs more human-friendly.