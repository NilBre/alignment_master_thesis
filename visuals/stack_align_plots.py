import re
import numpy as np
import matplotlib.pyplot as plt
import json

names = ["test_low_lumi", 
	"test_normal_lumi"]#,
#	"modules_low_lumi", 
#	"modules_normal_lumi"]

fig1 = plt.figure()
fig1.suptitle("GoodLongTracks, chi2 / dofs vs. nTracks, 10k events")
ax1 = fig1.add_subplot(121)

for run in names:
  f = open('../output/stack_alignment/{}/AlignmentResults/parsedlog.json'.format(run))
  data = json.load(f)

  vals, dofs = [], []
  for j in data['total_chi2_vals']:
    vals.append(j)
  for jj in data['total_chi2_dofs']:
    dofs.append(jj)

  vals, dofs = np.array(vals), np.array(dofs)
  chi2_per_dofs = vals / dofs

  T1, T2, T3 = [], [], []
  for i in data['FT/T1']['nTracks']:
    T1.append(float(i))
  for ii in data['FT/T2']['nTracks']:
    T2.append(float(ii))
  for iii in data['FT/T3']['nTracks']:
    T3.append(float(iii))

  nTracks = np.array(T1) + np.array(T2) + np.array(T3)

  zipped_list = zip(nTracks, chi2_per_dofs)
  sorted_pairs = sorted(zipped_list)
  tuples = zip(*sorted_pairs)
  tracks, chi2dofs = [list(tuple) for tuple in tuples]

  for i in range(len(tracks)):
    tracks[i] /= 1000

  ax1.plot(tracks, chi2dofs, label="{}".format(run))
  ax1.legend()
  ax1.grid(True)
  ax1.set_xlabel("nTracks x 1000")
  ax1.set_ylabel("chi2 / dofs")
  f.close()
fig1.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_{}.pdf'.format(run))

ax2 = fig1.add_subplot(122)
for run in names:
  f = open('../output/stack_alignment/{}/AlignmentResults/parsedlog.json'.format(run))
  data = json.load(f)

  vals, dofs = [], []
  for j in data['total_chi2_vals']:
    vals.append(j)
  for jj in data['total_chi2_dofs']:
    dofs.append(jj)

  vals, dofs = np.array(vals), np.array(dofs)
  chi2_per_dofs = vals / dofs

  T1, T2, T3 = [], [], []
  for i in data['FT/T1']['nTracks']:
    T1.append(float(i))
  for ii in data['FT/T2']['nTracks']:
    T2.append(float(ii))
  for iii in data['FT/T3']['nTracks']:
    T3.append(float(iii))

  nTracks = np.array(T1) + np.array(T2) + np.array(T3)

  for i in range(len(tracks)):
    tracks[i] /= 1000

  ax2.plot(np.linspace(0,9,10), chi2_per_dofs, label="{}".format(run))
  ax2.legend()
  ax2.grid(True)
  ax2.set_xlabel("iteration num")
  ax2.set_ylabel("chi2 / dofs")

  f.close()
fig1.savefig('outfiles/chi2_outfiles/chi2_vs_iter_{}.pdf'.format(run))
plt.show()
plt.clf()
##################################################

f1 = open('../output/stack_alignment/test_low_lumi/AlignmentResults/parsedlog.json')
f2 = open('../output/stack_alignment/test_normal_lumi/AlignmentResults/parsedlog.json')

data1 = json.load(f1)

fig2 = plt.figure()
fig2.suptitle("GoodLongTracks, nu=3.8 vs nu=7.6, 10k events")
ax1 = fig2.add_subplot(121)

vals, dofs = [], []
for j in data1['total_chi2_vals']:
  vals.append(j)
for jj in data1['total_chi2_dofs']:
  dofs.append(jj)

vals, dofs = np.array(vals), np.array(dofs)
chi2_per_dofs = vals / dofs

T1, T2, T3 = [], [], []
for i in data1['FT/T1']['nTracks']:
  T1.append(float(i))
for ii in data1['FT/T2']['nTracks']:
  T2.append(float(ii))
for iii in data1['FT/T3']['nTracks']:
  T3.append(float(iii))

nTracks = np.array(T1) + np.array(T2) + np.array(T3)

zipped_list = zip(nTracks, chi2_per_dofs)
sorted_pairs = sorted(zipped_list)
tuples = zip(*sorted_pairs)
tracks, chi2dofs = [list(tuple) for tuple in tuples]

for i in range(len(tracks)):
  tracks[i] /= 1000

ax1.plot(np.linspace(0,9,10), chi2dofs, label="low lumi, stationslayers")
ax1.legend()
ax1.grid(True)
ax1.set_xlabel("nTracks x 1000")
ax1.set_ylabel("chi2 / dofs")

f1.close()
###
ax2 = fig2.add_subplot(122)
data2 = json.load(f2)

vals, dofs = [], []
for j in data2['total_chi2_vals']:
  vals.append(j)
for jj in data2['total_chi2_dofs']:
  dofs.append(jj)

vals, dofs = np.array(vals), np.array(dofs)
chi2_per_dofs = vals / dofs

T1, T2, T3 = [], [], []
for i in data2['FT/T1']['nTracks']:
  T1.append(float(i))
for ii in data2['FT/T2']['nTracks']:
  T2.append(float(ii))
for iii in data2['FT/T3']['nTracks']:
  T3.append(float(iii))

nTracks = np.array(T1) + np.array(T2) + np.array(T3)

zipped_list = zip(nTracks, chi2_per_dofs)
sorted_pairs = sorted(zipped_list)
tuples = zip(*sorted_pairs)
tracks, chi2dofs = [list(tuple) for tuple in tuples]

for i in range(len(tracks)):
  tracks[i] /= 1000

ax2.plot(np.linspace(0, 9, 10), chi2dofs, label="normal lumi, stationslayers")
ax2.legend()
ax2.grid(True)
ax2.set_xlabel("nTracks x 1000")
ax2.set_ylabel("chi2 / dofs")

f2.close()
fig2.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_{}.pdf'.format(run))
plt.show()
plt.clf()

###################################
f3 = open('../output/stack_alignment/test_low_lumi/AlignmentResults/parsedlog.json')
f4 = open('../output/stack_alignment/test_normal_lumi/AlignmentResults/parsedlog.json')

data3 = json.load(f3)

fig3 = plt.figure()
fig3.suptitle("GoodLongTracks, nu=3.8 vs nu=7.6, 10k events")
ax1 = fig3.add_subplot(121)

vals, dofs = [], []
for j in data3['total_chi2_vals']:
  vals.append(j)
for jj in data3['total_chi2_dofs']:
  dofs.append(jj)

vals, dofs = np.array(vals), np.array(dofs)
chi2_per_dofs = vals / dofs

T1, T2, T3 = [], [], []
for i in data3['FT/T1']['nTracks']:
  T1.append(float(i))
for ii in data3['FT/T2']['nTracks']:
  T2.append(float(ii))
for iii in data3['FT/T3']['nTracks']:
  T3.append(float(iii))

nTracks = np.array(T1) + np.array(T2) + np.array(T3)

for i in range(len(tracks)):
  tracks[i] /= 1000

ax1.plot(np.linspace(0,9,10), chi2_per_dofs, label="low lumi, stationslayers")
ax1.legend()
ax1.grid(True)
ax1.set_xlabel("iteration num")
ax1.set_ylabel("chi2 / dofs")

f3.close()

#
ax2 = fig3.add_subplot(122)
data4 = json.load(f4)

vals, dofs = [], []
for j in data4['total_chi2_vals']:
  vals.append(j)
for jj in data4['total_chi2_dofs']:
  dofs.append(jj)

vals, dofs = np.array(vals), np.array(dofs)
chi2_per_dofs = vals / dofs

T1, T2, T3 = [], [], []
for i in data4['FT/T1']['nTracks']:
  T1.append(float(i))
for ii in data4['FT/T2']['nTracks']:
  T2.append(float(ii))
for iii in data4['FT/T3']['nTracks']:
  T3.append(float(iii))

nTracks = np.array(T1) + np.array(T2) + np.array(T3)

for i in range(len(tracks)):
  tracks[i] /= 1000

ax2.plot(np.linspace(0,9,10), chi2_per_dofs, label="normal lumi, stationslayers")
ax2.legend()
ax2.grid(True)
ax2.set_xlabel("iteration num")
ax2.set_ylabel("chi2 / dofs")

f4.close()

fig3.savefig('outfiles/chi2_outfiles/chi2_vs_iter_{}.pdf'.format(run))

plt.show()
plt.clf()
