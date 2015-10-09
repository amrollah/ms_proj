
from pysolar import radiation
from pysolar.solar import get_altitude, get_azimuth

latitude_deg = 42.3 # positive in the northern hemisphere
longitude_deg = -71.4 # negative reckoning west from prime meridian in Greenwich, England
d = datetime.datetime(2007, 2, 18, 15, 13, 1, 130320)
altitude_deg = get_altitude(latitude_deg, longitude_deg, d)
azimuth_deg = get_azimuth(latitude_deg, longitude_deg, d)
print(radiation.get_radiation_direct(d, altitude_deg))