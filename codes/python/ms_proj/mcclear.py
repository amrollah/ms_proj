from time import strptime, strftime
import csv

filename_mcclear = 'C:\data\mcclear\irradiation-clear_sky_cav_2015_till_september.csv'
mcclear_out = 'C:\data\mcclear\irradiation-clear_sky_cav_2015_till_september_old_v2.csv'

# t = strptime("2015-08-09T09:40:00.0", "%Y-%m-%dT%H:%M:%S.0")
# df = strftime("%d-%m-%Y %H:%M:%S", t)

with open(filename_mcclear, 'r') as mcfile:
    mc_reader = csv.reader(mcfile, delimiter=';')
    with open(mcclear_out, 'w', newline="") as outfile:
        mc_writer = csv.writer(outfile, delimiter=';')
        for row in mc_reader:
            t = strftime("%d-%m-%Y %H:%M:%S",strptime(row[0].split('/')[0], "%Y-%m-%dT%H:%M:%S.0")).split(' ')
            print(t)
            mc_writer.writerow(t + row[1:])

