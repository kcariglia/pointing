"""
docs here

K Cariglia, 2/26/25
"""

import argparse
import datetime
import json
import os, os.path

import astropy.coordinates
import astropy.time
import astropy.units
import numpy
import pandas

# westford site constants
siteLat = 42.612949 # degrees
siteLon = -71.493797 # degrees
siteEl = 86.764481 # meters

# assume satellite list + src_ods.json are in this directory
satListCsv = "/Users/cariglia/pointing/Starlink_DTC_culmination_Westford_2025-02-26T050000.000_2025-02-27T050000.000.csv"
srcJson = "/Users/cariglia/pointing/vo5057wf_ods.json"
satListDf = pandas.read_csv(satListCsv)
srcDf = pandas.read_json(srcJson)

# build list of satellite timestamps for easy parsing later
# indices correspond to those in satListDf
sat_timestamps = numpy.asarray([datetime.datetime.strptime(i, "%Y-%m-%dT%H:%M:%S.%f").timestamp() for i in satListDf['timestamp']])

# keep track of which satellites are near which source
# key = src_id, value = sat_list
src_sat_dict = {}
for s in range(len(srcDf)):
    # get needed src info
    src = srcDf['ods_data'][s]
    src_name = src['src_id']
    src_start = datetime.datetime.strptime(src['src_start_utc'], "%Y-%m-%dT%H:%M:%S.%f")
    src_end = datetime.datetime.strptime(src['src_end_utc'], "%Y-%m-%dT%H:%M:%S.%f")
    src_ra = src['src_ra_j2000_deg']
    src_dec = src['src_dec_j2000_deg']
    sat_list = []

    # get src az, el for both start and end times
    p = astropy.coordinates.SkyCoord(ra=src_ra, dec=src_dec, unit="deg")
    l = astropy.coordinates.EarthLocation(lat=siteLat, lon=siteLon, height=siteEl)
    t1 = astropy.time.Time(src_start)
    t2 = astropy.time.Time(src_end)
    azel1 = p.transform_to(astropy.coordinates.AltAz(obstime=t1, location=l))
    azel2 = p.transform_to(astropy.coordinates.AltAz(obstime=t2, location=l))

    # get indices where satellite overhead timestamp is within this src's start and end time
    potential_sat_idxs = numpy.argwhere((sat_timestamps >= src_start.timestamp()) & (sat_timestamps < src_end.timestamp()))

    for sat_idx in list(potential_sat_idxs):
        [sat_idx] = sat_idx
        sat_idx = int(sat_idx)
        sat = satListDf.iloc[sat_idx]
        sat_name = sat['sat']
        sat_el = sat['elevations']
        sat_az = sat['azimuths']

        print(f"(ut){datetime.datetime.fromtimestamp(sat_timestamps[sat_idx])}  (sel){sat_el}  (el1){azel1.alt.value} (el2){azel2.alt.value} (saz){sat_az}  (az1){azel1.az.value}  (az2){azel2.az.value}")

        f1 = abs(azel1.az.value - sat_az) <= 10
        f2 = abs(azel2.az.value - sat_az) <= 10
        f3 = abs(azel1.alt.value - sat_el) <= 10
        f4 = abs(azel2.alt.value - sat_el) <= 10
        if (f1 or f2):
            if (f3 or f4):
                # this satellite is within 10 degrees of az, el to the source
                # add to sat_list
                sat_list.append(sat_name)

    # remove duplicates and add to src_sat_dict
    sat_list = list(set(sat_list))
    src_sat_dict[src_name] = sat_list


with open("src_sat_list.txt", "w") as f:
    f.write("# src_id       sat_list\n")
    for k in src_sat_dict.keys():
        writestr = k + "      " + ",".join(src_sat_dict[k]) + "\n"



    