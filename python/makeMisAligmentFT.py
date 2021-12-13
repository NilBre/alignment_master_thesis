__doc__ = ''' File for submit misalign the FT detector '''


import argparse
from shutil import copyfile
import os


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('level', choices=["Layer", "Module", "Mats"],
                    help='Which level you want misalign?')
parser.add_argument('type', choices=["Translation", "Rotation", "all"],
                    help='Which type of tranformation?')
parser.add_argument('-t', type=float, default=10,
                    help='Misaligment for translations in microns (default: 10um)')
parser.add_argument('-r', type=float, default=10,
                    help='Misaligment for rotations in millradiant (default: 10mrad)')
parser.add_argument('--taxis', choices=['x', 'y', 'z', 'xz', 'all'], default='all',
                    help='Which axis to translate (default: all)')
parser.add_argument('--raxis', choices=['x', 'y', 'z', 'all'], default='all',
                    help='Which axis to rotate (default: all)')
parser.add_argument('--reset', action='store_true',
                    help='generate a xml with all zeros')
args = parser.parse_args()


#DBDIR = os.getenv('GITCONDDBPATH')
DBDIR = "./nominal_xml"
dlevel = {
     "Layer" : 'FTSystem.xml'
    ,"Module": 'Modules.xml'
    ,"Mats"  : 'Mats.xml'
}
#myfile = DBDIR + '/SIMCOND/Conditions/FT/Alignment/%s'%dlevel[args.level]
myfile = DBDIR + '/%s'%dlevel[args.level]
copyfile(myfile, dlevel[args.level])
print('Input file:', myfile)



#create Layer strings
stations = ["T1", "T2", "T3"]
layers = ["X1", "U", "V", "X2"]
#layers = ["X1"]
quarters = ['Q0', 'Q1', 'Q2', 'Q3']
modules = ['M0', 'M1', 'M2', 'M3', 'M4']
glayers = []
for s in stations:
  for l in layers:
    for q in quarters:
      for m in modules:
        glayers.append(s + l + q + m)
print('Misalign only modules/mats in Layer\n', glayers)

import numpy as np
rnd = np.random
#-------------------------------------------------------------------------------
mu, sigma = 0., args.t*1e-3 #mm
mu_r, sigma_r = 0., args.r*1e-3 #rad
axis_to_idx = {
   'x': [1,0,0], 
   'y': [0,1,0], 
   'z': [0,0,1], 
   'xz': [1,0,1], 
   'all': [1,1,1],
   }
mask_axis_t = np.asarray( axis_to_idx[args.taxis] )
mask_axis_r = np.asarray( axis_to_idx[args.raxis] )
if args.reset: mask_axis = np.asarray( [0,0,0] )
#-------------------------------------------------------------------------------


import xml.etree.ElementTree as ET

tree = ET.parse(dlevel[args.level])
root = tree.getroot()

for cond in root.iter('condition'):
    #print(cond.attrib)
    skip_this = True
    for glay in glayers:
      if glay in cond.attrib['name']:
        skip_this = False

    if skip_this:
        continue


    for pvec in cond.iter('paramVector'):
      if pvec.attrib['name'] == 'dPosXYZ' and ('Translation' == args.type or 'all' == args.type):
        random_shift = rnd.normal(mu, sigma, 3)*mask_axis_t
        pvec.text = "{0} {1} {2}".format(random_shift[0], random_shift[1], random_shift[2])
      if pvec.attrib['name'] == 'dRotXYZ' and ('Rotation' == args.type or 'all' == args.type):
        random_shift = rnd.normal(mu_r, sigma_r, 3)*mask_axis_r
        pvec.text = "{0} {1} {2}".format(random_shift[0], random_shift[1], random_shift[2])

tree.write(dlevel[args.level])

with open(dlevel[args.level], "r+") as f:
  old = f.read() # read everything in the file
  f.seek(0) # rewind
  f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n<!DOCTYPE DDDB SYSTEM "git:/DTD/structure.dtd">\n' + old) # write the new line before


print('Output file: ', dlevel[args.level])

#Add this at the begin of the file
#<?xml version="1.0" encoding="ISO-8859-1"?>
#<!DOCTYPE DDDB SYSTEM "git:/DTD/structure.dtd">
