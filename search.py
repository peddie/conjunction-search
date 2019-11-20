#!/usr/bin/env python3

import numpy as np
import ephem  # for magnitude calculations
from datetime import datetime, timedelta, date
from skyfield.api import Star, load, Topos, utc
from skyfield.data import hipparcos
from skyfield.units import Angle
from skyfield import almanac
import pytz
from tzwhere import tzwhere
import re
import time
from tqdm import tqdm
import argparse


def datetime_range(start, end, delta):
    current = start
    while current < end:
        yield current
        current += delta


def ll_string_to_float(string):
    cardinal = ''.join(c for c in string if c in "nsewNSEW")
    magnitude = ''.join(c for c in string if c.isdigit() or c == '.')
    sign = 1
    if cardinal in "sSwW":
        sign = -1

    return sign * float(magnitude)


def float_to_lat_string(value):
    cardinal = 'N'
    if value < 0:
        cardinal = 'S'

    magnitude = str(abs(value))
    return f"{magnitude} {cardinal}"


def float_to_lon_string(value):
    cardinal = 'E'
    if value < 0:
        cardinal = 'W'

    magnitude = str(abs(value))
    return f"{magnitude} {cardinal}"


def parse_date(s):
    return datetime.strptime(s, '%Y-%m-%d')


def angular_distance(va, vb):
    dlambda = np.fabs(vb[1].radians - va[1].radians)
    phi1 = va[0].radians
    phi2 = vb[0].radians

    s1 = np.sin(phi1)
    s2 = np.sin(phi2)

    c1 = np.cos(phi1)
    c2 = np.cos(phi2)

    slam = np.sin(dlambda)
    clam = np.cos(dlambda)

    return np.arctan(np.sqrt((c2*slam)**2 + (c1*s2 - s1*c2*clam)**2) / (s1*s2 + c1*c2*clam))


def load_stars(max_magnitude):
    """Load the Hipparcos catalog data and find all applicable stars.  Returns a
    list of StarObs, one per star.
    """
    print("[  CATALOG] Loading star catalog data.")
    with load.open(hipparcos.URL) as f:
        df = hipparcos.load_dataframe(f)

        df = df[df['magnitude'] <= max_magnitude]
        df = df[df['ra_degrees'].notnull()]

        print(f'[  CATALOG] Found {len(df)} stars brighter than magnitude {max_magnitude:.4f}.')
        return [StarObs(str(idx), Star.from_dataframe(df.loc[idx]), df.loc[idx].magnitude)
                for idx in df.index]


class Observation:
    def __init__(self, time, site, obj, magnitude, radius):
        self.t = time
        self.site = site
        self.obj = obj
        # Compute right at the get-go
        self.radius = radius
        self.magnitude = magnitude
        # This is also precomputed due to the location precomputation thing
        self.apparent_alt_az = self.site.apparent_alt_az(self.t, self.obj)

    def alt_az(self):
        return self.apparent_alt_az

    def apparent_magnitude(self):
        return self.magnitude

    def apparent_radius(self):
        return self.radius

    def time(self):
        return self.t

    def summarize(self):
        altaz = self.alt_az()
        magnitude = self.apparent_magnitude()
        radius = self.apparent_radius()

        visibility = 'visible' if self.visible() else 'not visible'

        return f"""{self.obj.name}, at {self.t} {self.site.timezone}:
  altitude {altaz[0]}
  azimuth  {altaz[1]}
  magnitude {magnitude}
  radius {radius} radians
  {visibility} from {self.site}
"""

    def visible(self):
        return self.alt_az()[0].radians > 0

    def body(self):
        return self.obj


class LocationCacheError(Exception):
    def __init__(self, cached_time, current_time):
        self.cached_time = cached_time
        self.current_time = current_time

    def __str__(self):
        return f"""LocationCacheError:
Cached location is for {self.cached_time}
Attempted observation is for {self.current_time}
Please call 'set_t()' to update the location cache before observing!
"""
    def __repr__(self):
        return self.__str__()


class Site:
    def __init__(self, lat, lon, name):
        self.topos = Topos(lat, lon)
        self.lat = lat
        self.lon = lon
        self.ephemeris = load('de430.bsp')
        self.location = self.ephemeris['earth'] + self.topos
        self.ts = load.timescale()
        self.name = name
        lat_value, lon_value = ll_string_to_float(lat), ll_string_to_float(lon)
        self.timezone = pytz.timezone(tzwhere.tzwhere().tzNameAt(lat_value, lon_value))
        print(f"Time zone at {name}: {self.timezone}")

    def localize(self, some_time):
        return self.timezone.localize(some_time)

    def set_t(self, t):
        """Hack to reduce computation time.  Whenever we start computing conjunctions
        for a new time t, we precompute our location here.  Then in
        `apparent_alt_az()`, we use this cached location object.  This is
        extremely ugly, but it results in a >3x speedup in conjunctions per
        second.

        """
        self.location_at_t = self.location.at(self.ts.utc(t))
        self.cached_t = t

    def apparent_alt_az(self, t, obj):
        """See the hack documented in `set_t()` to learn why `t` is not used here.

        """
        if t is not self.cached_t:
            raise LocationCacheError(self.cached_t, t)

        alt, az, _ = self.location_at_t.observe(obj.skyfield_obj).apparent().altaz()
        return (alt, az)

    def observe(self, t, obj):
        obj.update_state(t)
        return Observation(t, self, obj, obj.get_magnitude(), obj.get_radius())

    def next_set_rise(self, day):
        # By adding 2 days, we can be sure to find today's sunset and the following (tomorrow
        # morning's) sunrise.
        #
        # TODO(MP): Like all other times, the UTC situation here is a bit crappy.  We should really
        # just start searching at noon local time.
        begin = self.ts.utc(day)
        end = self.ts.utc(day + timedelta(days=2))
        times, ups = almanac.find_discrete(begin, end,
                                           almanac.sunrise_sunset(self.ephemeris, self.topos))
        times = [t.astimezone(self.timezone) for t in times]

        # Find the first sunset
        sunset = next(horizon for horizon in zip(times, ups) if not horizon[1])[0]
        # Find the first sunrise after the first sunset
        sunrise = next(horizon for horizon in zip(times, ups) if horizon[1] and horizon[0] > sunset)[0]
        return (sunset, sunrise)

    def __str__(self):
        return f'{self.name}, at {self.lat}, {self.lon} (Earth)'

    def __repr__(self):
        return self.__str__()


class SolarSystemBody():
    def __init__(self, name, skyfield_obj, ephem_obj):
        self.name = name
        self.skyfield_obj = skyfield_obj
        self.ephem_obj = ephem_obj

    def __str__(self):
        return f'SolarSystemBody({self.name})'

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return self.__str__() <= other.__str__()

    def update_state(self, t):
        self.ephem_obj.compute(t)
        return None

    def get_magnitude(self):
        return self.ephem_obj.mag

    def get_radius(self):
        return self.ephem_obj.radius


def load_solar_system():
    print("[EPHEMERIS] Loading ephemeris data for planets / moon / Pluto.")
    planet_barycenters = load('de430.bsp')

    return [
        SolarSystemBody('Mercury', planet_barycenters['mercury'], ephem.Mercury()),
        SolarSystemBody('Venus', planet_barycenters['venus'], ephem.Venus()),
        SolarSystemBody('Moon', planet_barycenters['moon'], ephem.Moon()),
        SolarSystemBody('Mars', load('mar097.bsp')['mars'], ephem.Mars()),
        SolarSystemBody('Jupiter', load('jup310.bsp')['jupiter'], ephem.Jupiter()),
        SolarSystemBody('Saturn', load('sat319.bsp')['saturn'], ephem.Saturn()),
        SolarSystemBody('Uranus', load('ura111.bsp')['uranus'], ephem.Uranus()),
        SolarSystemBody('Neptune', load('nep081.bsp')['neptune'], ephem.Neptune()),
        SolarSystemBody('Pluto', load('plu055.bsp')['pluto'], ephem.Pluto()),
    ]


class StarObs():
    def __init__(self, name, skyfield_obj, magnitude):
        self.name = name
        self.skyfield_obj = skyfield_obj
        self.magnitude = magnitude

    def __str__(self):
        return f'StarObs({self.name}, mag {self.magnitude})'

    def __repr__(self):
        return self.__str__()

    def update_state(self, t):
        return None

    def get_magnitude(self):
        return self.magnitude

    def get_radius(self):
        return 0.0


class Constraint():
    def __init__(self, max_magnitude, max_distance, max_magnitude_delta):
        self.max_magnitude = max_magnitude
        self.max_distance = max_distance
        self.max_magnitude_delta = max_magnitude_delta


class Conjunction():
    def __init__(self, first_obs, second_obs, constraint):
        self.first_obs = first_obs
        self.second_obs = second_obs
        assert first_obs.time() == second_obs.time()
        self.time = first_obs.time()
        self.constraint = constraint

    def valid(self):
        if not self.first_obs.visible() or not self.second_obs.visible():
            return False

        maga, magb = self.first_obs.apparent_magnitude(), self.second_obs.apparent_magnitude()
        if not maga <= self.constraint.max_magnitude or not magb <= self.constraint.max_magnitude:
            return False

        if np.abs(maga - magb) > self.constraint.max_magnitude_delta:
            return False

        if np.abs(self.__distance()) > self.constraint.max_distance:
            return False

        return True


    def __distance(self):
        # In the unlikely event of overlapping objects, we'll get a negative distance.
        return (angular_distance(self.first_obs.alt_az(),
                                 self.second_obs.alt_az()) -
                self.first_obs.apparent_radius() -
                self.second_obs.apparent_radius())

    def __magnitude_difference(self):
        first_mag, second_mag = self.first_obs.apparent_magnitude(), self.second_obs.apparent_magnitude()
        return np.abs(first_mag - second_mag)

    def __max_magnitude(self):
        first_mag, second_mag = self.first_obs.apparent_magnitude(), self.second_obs.apparent_magnitude()
        return max(first_mag, second_mag)

    def score(self):
        distance = self.__distance()
        magnitude_difference = self.__magnitude_difference()
        max_magnitude = self.__max_magnitude()

        return distance + magnitude_difference + max_magnitude

    def __str__(self):
        return f"""Conjunction between {self.first_obs.obj.name} and {self.second_obs.obj.name}:
                time: {self.time}
            distance: {self.__distance() * 180 / np.pi:.3f} degrees
magnitude difference: {self.__magnitude_difference():.3f} ({np.power(2.152, self.__magnitude_difference()):.3f} times brighter)
   maximum magnitude: {self.__max_magnitude():.3f}
               score: {self.score():.5f}
"""

    def __repr__(self):
        return self.__str__()

    def bodies(self):
        # Return bodies in normalized order
        s = sorted([self.first_obs.body(), self.second_obs.body()])
        return (s[0], s[1])


def sort_conjunctions_by_score(conjunctions):
    return sorted(conjunctions, key=lambda x: x.score())


def scan_night(site, day, constraint):
    t0 = time.perf_counter()
    dt = timedelta(minutes=30)
    sunset, sunrise = site.next_set_rise(day)

    conjunctions = {}
    num_conjunctions = len(solar_system) * len(sky_objects)
    for t in tqdm(datetime_range(sunset, sunrise, dt),
                  desc='Time windows', unit='win', leave=None,
                  total=len(list(datetime_range(sunset, sunrise, dt)))):
        site.set_t(t)
        # Form all the possible conjunction pairs
        conj = [
            c for c in tqdm((Conjunction(site.observe(t, p), site.observe(t, q), constraint)
                             # Scan only solar system objects for the first index, since stars
                             # are obviously static
                             for p in solar_system
                             for q in sky_objects
                             if p != q),
                            desc='Conjunctions scanned', unit='conj', leave=None,
                            total=num_conjunctions)
            if c.valid()

        ]

        # Filter out invalid conjunctions
        # conj = [c for c in conj if c.valid()]
        for c in conj:
            # Save to the master dictionary, indexed by the bodies
            if c.bodies() in conjunctions:
                conjunctions[c.bodies()].append(c)
            else:
                conjunctions[c.bodies()] = [c]

    outputs = []
    for _, conjs in conjunctions.items():
        # Find the best output for the whole night for this pair of bodies
        outputs.append(sort_conjunctions_by_score(conjs)[0])

    t1 = time.perf_counter()
    return sort_conjunctions_by_score(outputs), t1 - t0


def scan_days(location, t0_, t1_, constraint=Constraint(5.0, 1 * np.pi / 180, 3.0)):
    t0 = location.localize(t0_)
    t1 = location.localize(t1_)

    dday = timedelta(days=1)

    best_conjunctions = {}

    print(f"[   SEARCH] Scanning nights from {t0} to {t1}")
    print(f"[   SEARCH] This will take a while!")
    calc_times = []
    for day in tqdm(datetime_range(t0, t1, dday),
                    desc='Nights', unit='night', leave=None,
                    total=len(list(datetime_range(t0, t1, dday)))):
        best_conjunctions[day], calc_time = scan_night(location, day, constraint)
        calc_times.append(calc_time)

    mean_time_per_night = np.mean(calc_times)
    total_time = np.sum(calc_times)
    print(f"[   SEARCH] Scanning took {total_time:.3f}s for {len(calc_times)} nights.")
    print(f"[   SEARCH] Mean scan time was {mean_time_per_night:.3f}s per night.")

    return best_conjunctions


def display_results(results, site, t0, t1, best=5):
    print(f"""Displaying events in {site} between
    {t0}
 and
    {t1}
 :
""")
    for night, conjunctions in results.items():
        num_conjunctions = len(conjunctions)
        if num_conjunctions == 0:
            continue
        print(night, f'({num_conjunctions} event{"" if num_conjunctions == 1 else "s"})')
        if num_conjunctions < best:
            print(sort_conjunctions_by_score(conjunctions))
        else:
            print(sort_conjunctions_by_score(conjunctions)[0:best])


def get_args():
    parser = argparse.ArgumentParser(description=f"""
Search for astronomical conjunctions.

By default:
 - location is Brisbane, QLD, Australia.
 - start date is tonight ({date.today()}).
 - end date is one week from now ({date.today() + timedelta(days=7)}).

Other default values are displayed below in parentheses.
""", formatter_class=argparse.RawTextHelpFormatter)
    ## Observation setup
    obs_options = parser.add_argument_group('OBSERVER CONFIGURATION')
    obs_options.add_argument('-a', '--latitude', metavar='LAT', default='27.49665 S',
                             help='Latitude of the observation location ("%(default)s")')
    obs_options.add_argument('-o', '--longitude', metavar='LON', default='152.9883 E',
                             help='Longitude of the observation location ("%(default)s")')
    obs_options.add_argument('--site-name', metavar='NAME',
                             default='Home',
                             help='Name of the observing site (just makes output prettier) (%(default)s)')

    ## Search setup
    search_options = parser.add_argument_group('SEARCH CONFIGURATION')
    search_options.add_argument('--start-date', metavar='START',
                             type=parse_date, default=date.today(),
                             help='Date of the first night (sunset) to search (%(default)s)')
    group = search_options.add_mutually_exclusive_group()
    group.add_argument('--num-nights', metavar='NIGHTS',
                       type=int, default=7,
                       help='Number of nights after the start date to search (%(default)s)')
    group.add_argument('--end-date', metavar='END',
                       type=parse_date,
                       help='Date of the last sunrise up to which to search (calculated from --num-nights)')

    search_options.add_argument('--max-star-magnitude', metavar='MAXSTARMAG',
                                type=float, default=4.0,
                                help='Maximum star magnitude to consider when searching (%(default)s)')
    search_options.add_argument('--max-magnitude', metavar='MAXMAG',
                                type=float, default=5.0,
                                help='Maximum planet/moon magnitude to consider when searching (%(default)s)')
    search_options.add_argument('--max-magnitude-delta', metavar='DELTA',
                                type=float, default=3.0,
                                help='Maximum magnitude discrepancy (between objects in a conjunction) to consider when searching (%(default)s)')
    search_options.add_argument('--max-angle', metavar='MAXANGLE',
                                type=float, default=1.0,
                                help='Maximum planet/moon magnitude to consider when searching in degrees (%(default)s)')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()

    solar_system = load_solar_system()
    stars = load_stars(max_magnitude=args.max_star_magnitude)

    sky_objects = solar_system + stars

    brisbane = Site(args.latitude, args.longitude, args.site_name)

    t0 = datetime.combine(args.start_date, datetime.min.time())
    t1 = datetime.combine(args.end_date
                          if args.end_date
                          else args.start_date + timedelta(days=args.num_nights),
                          datetime.min.time())

    results = scan_days(brisbane, t0, t1,
                        constraint=Constraint(args.max_magnitude,
                                              args.max_angle * np.pi / 180,
                                              args.max_magnitude_delta))
    display_results(results, brisbane, t0, t1)
