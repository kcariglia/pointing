"""
---------------------------
Parses a Westford telescope pointing file and outputs it in the ODS JSON format as described in the link below:
https://www.seti.org/hcro/ods
---------------------------
---------------------------
Assumes there exists a local copy of my fork of David Doeber's "odsutils" library that is installed in the 
active python environment. If not, please run the following:
> git clone https://github.com/kcariglia/odsutils.git
and install in the python environment:
> cd odsutils
> pip install .
---------------------------
Also assumes the following are installed in the active python environment:
-- python 3.12
-- astropy
---------------------------
---------------------------
Usage:
python parse_pointing.py input_pointing_file.snp
---------------------------
---------------------------
K Cariglia, 1/3/25
"""

import argparse
import datetime
import json
import os, os.path

import astropy.coordinates
import astropy.time
import astropy.units

# hardcoded path to github checkout of odsutils/scripts
ODS_PTH = "/Users/cariglia/odsutils/scripts"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="pointing file to parse", default=None)
    args = vars(parser.parse_args())

    if not args["filename"]:
        raise IOError("No pointing file given")
    
    fname = args["filename"]

    with open(fname, "r") as f:
        lines = f.readlines()

    # assume these are constant
    siteID = "Wf"
    siteLat = 42.612949 # degrees
    siteLon = -71.493797 # degrees
    siteEl = 86.764481 # meters
    radius = .3 # degrees
    slew_az = 200.0 / 60.0 # degrees per second, az
    slew_el = 2.0 # degrees per second, el

    # need to reset these when there is a new source
    src_name = None
    src_is_pulsar = False
    cor_integ_time = 1
    src_ra = None
    src_dec = None

    startTime = None
    lastTime = None
    slewTime = None
    recordLen = None
    trk_rate_dec = 0
    trk_rate_ra = 0
    freq_lower = 2000000000
    freq_upper = 14000000000
    notes = ''

    # store all needed info in a list of dicts
    source_list = []
    this_source = {}

    # separate header info from the multiple sources just in case
    for line in lines:
        if line.startswith("\""):
            # part of the header, dont need any info here
            pass
        elif line.startswith("!"):
            # timestamp 
            lastTime = datetime.datetime.strptime(line[1:], "%Y.%j.%H:%M:%S\n")
        elif line.startswith("source="):
            # denotes new source

            if lastTime is not None:
                # need to figure out start time
                # assume we have last time

                if len(source_list) > 0:
                    # we need to calculate slew time

                    # point 1 is prev point ra, dec, prev point lastTime
                    p1 = astropy.coordinates.SkyCoord(ra=source_list[-1]['src_ra_j2000_deg'], dec=source_list[-1]['src_dec_j2000_deg'], unit="deg")
                    l1 = astropy.coordinates.EarthLocation(lat=siteLat, lon=siteLon, height=siteEl)
                    t1 = astropy.time.Time(source_list[-1]['src_end_utc'])

                    azel1 = p1.transform_to(astropy.coordinates.AltAz(obstime=t1, location=l1))

                    # point 2 is current point ra, dec, lastTime - (recordLen + 4)
                    p2 = astropy.coordinates.SkyCoord(ra=src_ra, dec=src_dec, unit="deg")
                    t2 = astropy.time.Time(lastTime - datetime.timedelta(seconds=(recordLen+4)))

                    azel2 = p2.transform_to(astropy.coordinates.AltAz(obstime=t2, location=l1))

                    dif_az = azel2.az.value - azel1.az.value
                    dif_el = azel2.alt.value - azel1.alt.value
                    az_time = abs(dif_az) / slew_az
                    el_time = abs(dif_el) / slew_el
                    slewTime = az_time + el_time

                # times should be contiguous
                if len(source_list) == 0:
                    startTime = lastTime - datetime.timedelta(seconds=(recordLen + 4 + slewTime))
                else:
                    startTime = datetime.datetime.strptime(source_list[-1]['src_end_utc'], "%Y-%m-%dT%H:%M:%S.%f")

                # populate dict for this source
                this_source['site_id'] = siteID
                this_source['site_lat_deg'] = siteLat
                this_source['site_lon_deg'] = siteLon
                this_source['site_el_m'] = siteEl
                this_source['src_id'] =  src_name
                this_source['src_is_pulsar_bool'] = src_is_pulsar
                this_source['corr_integ_time_sec'] = cor_integ_time
                this_source['src_ra_j2000_deg'] = src_ra
                this_source['src_dec_j2000_deg'] = src_dec
                this_source['src_radius'] = radius
                this_source['src_start_utc'] = startTime.strftime("%Y-%m-%dT%H:%M:%S.%f")
                this_source['src_end_utc'] = lastTime.strftime("%Y-%m-%dT%H:%M:%S.%f")
                this_source['slew_sec'] = slewTime
                this_source['trk_rate_dec_deg_per_sec'] = trk_rate_dec
                this_source['trk_rate_ra_deg_per_sec'] = trk_rate_ra
                this_source['freq_lower_hz'] = freq_lower
                this_source['freq_upper_hz'] = freq_upper
                this_source['notes'] = notes
                source_list.append(this_source)

                # reset fields for this source dict
                this_source = {}
                src_name = None
                src_is_pulsar = False
                cor_integ_time = 1
                src_ra = None
                src_dec = None

                startTime = None
                lastTime = None
                slewTime = None
                recordLen = None
                trk_rate_dec = 0
                trk_rate_ra = 0
                freq_lower = 2000000000
                freq_upper = 14000000000
                notes = ''

                # get source name, ra, dec
                src_name = line[7:].split(',')[0]
                src_ra_str = line[7:].split(',')[1]
                src_dec_str = line[7:].split(',')[2]

                # assume ra str is always has the same number of characters
                src_ra = astropy.coordinates.Angle(src_ra_str[0:2] + ':' + src_ra_str[2:4] + ':' + src_ra_str[4:], unit=astropy.units.hourangle).to_value(astropy.units.degree)
                
                if src_dec_str.startswith("-"):
                    src_dec = astropy.coordinates.Angle("{}d{}'{}\"".format(src_dec_str[0:3], src_dec_str[3:5], src_dec_str[5:])).to_value(astropy.units.degree)
                else:
                    src_dec = astropy.coordinates.Angle("{}d{}'{}\"".format(src_dec_str[0:2], src_dec_str[2:4], src_dec_str[4:])).to_value(astropy.units.degree)
                
            else:
                # this is the first time we see a source line, no slew time info yet
                slewTime = 0

                # get source name, ra, dec
                src_name = line[7:].split(',')[0]
                src_ra_str = line[7:].split(',')[1]
                src_dec_str = line[7:].split(',')[2]

                # assume ra str is always has the same number of characters
                src_ra = astropy.coordinates.Angle(src_ra_str[0:2] + ':' + src_ra_str[2:4] + ':' + src_ra_str[4:], unit=astropy.units.hourangle).to_value(astropy.units.degree)

                if src_dec_str.startswith("-"):
                    src_dec = astropy.coordinates.Angle("{}d{}'{}\"".format(src_dec_str[0:3], src_dec_str[3:5], src_dec_str[5:])).to_value(astropy.units.degree)
                else:
                    src_dec = astropy.coordinates.Angle("{}d{}'{}\"".format(src_dec_str[0:2], src_dec_str[2:4], src_dec_str[4:])).to_value(astropy.units.degree)


        elif line.startswith("mk6=record="):
            # denotes length of recording
            recordLen = int(line.split(":")[1])

        elif line.startswith("sched_end"):

            if len(source_list) > 0:
                    # we need to calculate slew time

                    # point 1 is prev point ra, dec, prev point lastTime
                    p1 = astropy.coordinates.SkyCoord(ra=source_list[-1]['src_ra_j2000_deg'], dec=source_list[-1]['src_dec_j2000_deg'], unit="deg")
                    l1 = astropy.coordinates.EarthLocation(lat=siteLat, lon=siteLon, height=siteEl)
                    t1 = astropy.time.Time(source_list[-1]['src_end_utc'])

                    azel1 = p1.transform_to(astropy.coordinates.AltAz(obstime=t1, location=l1))

                    # point 2 is current point ra, dec, lastTime - (recordLen + 4)
                    p2 = astropy.coordinates.SkyCoord(ra=src_ra, dec=src_dec, unit="deg")
                    t2 = astropy.time.Time(lastTime - datetime.timedelta(seconds=(recordLen+4)))

                    azel2 = p2.transform_to(astropy.coordinates.AltAz(obstime=t2, location=l1))

                    dif_az = azel2.az.value - azel1.az.value
                    dif_el = azel2.alt.value - azel1.alt.value
                    az_time = abs(dif_az) / slew_az
                    el_time = abs(dif_el) / slew_el
                    slewTime = az_time + el_time

            startTime = datetime.datetime.strptime(source_list[-1]['src_end_utc'], "%Y-%m-%dT%H:%M:%S.%f")

            # very last source of this file, add last dict to list
            this_source['site_id'] = siteID
            this_source['site_lat_deg'] = siteLat
            this_source['site_lon_deg'] = siteLon
            this_source['site_el_m'] = siteEl
            this_source['src_id'] =  src_name
            this_source['src_is_pulsar_bool'] = src_is_pulsar
            this_source['corr_integ_time_sec'] = cor_integ_time
            this_source['src_ra_j2000_deg'] = src_ra
            this_source['src_dec_j2000_deg'] = src_dec
            this_source['src_radius'] = radius
            this_source['src_start_utc'] = startTime.strftime("%Y-%m-%dT%H:%M:%S.%f")
            this_source['src_end_utc'] = lastTime.strftime("%Y-%m-%dT%H:%M:%S.%f")
            this_source['slew_sec'] = slewTime
            this_source['trk_rate_dec_deg_per_sec'] = trk_rate_dec
            this_source['trk_rate_ra_deg_per_sec'] = trk_rate_ra
            this_source['freq_lower_hz'] = freq_lower
            this_source['freq_upper_hz'] = freq_upper
            this_source['notes'] = notes
            source_list.append(this_source)

    # dump data into json file to set default values
    with open("ods_in.json", "w") as f:
        json.dump({"ods_data":source_list}, f)
    
    # ods_in.json -> odsuser.py -> ods_out.json
    outfile = args['filename'].split('.')[0] + "_ods.json"
    cmd = "$(which python) {} -o ods_in.json -d from_ods -w {}".format(os.path.join(ODS_PTH, "odsuser.py"), outfile)
    ret = os.system(cmd)
    if (ret != 0):
        raise Exception("Problem creating ODS file")

    # double check that we can still read this file with odsuser.py
    cmd = "$(which python) {} -o {} --view".format(os.path.join(ODS_PTH, "odsuser.py"), outfile)
    ret = os.system(cmd)
    if (ret != 0):
        raise Exception("Problem reading ODS file")

    # no problems, so get rid of ods_in.json
    if os.access("ods_in.json", os.R_OK) and os.access(outfile, os.R_OK):
        os.remove("ods_in.json")


        
