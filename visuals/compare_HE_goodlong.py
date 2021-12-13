import re
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.stats import sem

# data input
GLtracks = open("../../../output_alignment/20211206/goodlongtracks/config5_2019/AlignmentResults/parsedlog.json")
HEtracks = open("../../../output_alignment/20211206/highmomentumttracks/config5_2019/AlignmentResults/parsedlog.json")
f1 = open('../output/alignment_runs/done_runs/HE_runs/13/AlignmentResults/parsedlog.json')

GL_data = json.load(GLtracks)
HE_data = json.load(HEtracks)
num13_data = json.load(f1)

# total number of tracks = sum over stations and tracks
GL_T1, GL_T2, GL_T3 = [], [], []
HE_T1, HE_T2, HE_T3 = [], [], []
fT1, fT2, fT3 = [], [], []

stations = ["FT/T1", "FT/T2", "FT/T3"]

for s in stations:
  for i in GL_data[s]['nTracks']:
    if "1" in s:    
      GL_T1.append(float(i))
    elif "2" in s:
      GL_T2.append(float(i))
    else:
      GL_T3.append(float(i))
  for j in HE_data[s]['nTracks']:
    if "1" in s:
      HE_T1.append(float(j))
    elif "2" in s:
      HE_T2.append(float(j))
    else:
      HE_T3.append(float(j))
  for k in num13_data[s]['nTracks']:
    if "1" in s:
      fT1.append(float(k))
    elif "2" in s:
      fT2.append(float(k))
    else:
      fT3.append(float(k))

GL_nTracks = np.array(GL_T1) + np.array(GL_T2) + np.array(GL_T3)
HE_nTracks = np.array(HE_T1) + np.array(HE_T2) + np.array(HE_T3)

num13_nTracks = np.array(fT1) + np.array(fT2) + np.array(fT3)
# fehler des mittelerts: std / sqrt(n)
sem_GL = sem(GL_nTracks)
sem_HE = sem(HE_nTracks)

sem_num13 = sem(num13_nTracks)
# 3 digits
limit_GL = "{:.3f}".format(sem_GL)
limit_HE = "{:.3f}".format(sem_HE)
limit_num13 = "{:.3f}".format(sem_num13)

print(sem_GL, sem_HE, sem_num13)

print("GoodLongTracks: ", GL_nTracks)
print("HighMomentumTTracks: ", HE_nTracks)
print("HE, config5, 2020, 20 iterations: ", num13_nTracks)
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print("quotient HE / GL: ", GL_nTracks / HE_nTracks)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax1.errorbar(np.linspace(1, 10, 10), GL_nTracks, yerr=sem_GL, fmt="r-", label="GoodLongTracks: err={}".format(limit_GL))
ax1.legend()
ax1.grid(True)
ax1.set_xlabel('iteration')
ax1.set_ylabel('number of tracks')

ax2 = fig.add_subplot(122)
ax2.errorbar(np.linspace(1, 10, 10), HE_nTracks, yerr=sem_HE, fmt="b-", label="HighMomentumTTracks: err={}".format(limit_HE))
ax2.legend()
ax2.grid(True)
ax2.set_xlabel('iteration')
ax2.set_ylabel('number of tracks')
plt.show()
plt.clf()

plt.errorbar(np.linspace(1, 20, 20), num13_nTracks, yerr=sem_num13, fmt="g-", label="HE, c5, 2020, iter20: err={}".format(limit_num13))
plt.legend()
plt.grid()
plt.xlabel('iteration')
plt.ylabel('number of tracks')
plt.show()
