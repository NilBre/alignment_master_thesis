import re
import numpy as np
import matplotlib.pyplot as plt
import json

# loop over the 4 HE runs and without the 30k events
# names_30k = ["2019,TxTz,30k", "2019,TxTz,withModules,30k"]
names = ['2020, config5, TxTzRxRz', '2020, config5, TxRz', '2020, withModules, TxTzRxRz', '2020, withModules, TxRz']#, "2020, 16iter TxRz", "2019,TxRz,normal,6.5k", "2019,TxRz,modules,6.5k", "onlyModules,2019,TxRz", "onlyModules,2020,TxRz"]
names_2020 = ["config5, TxRz", "config5, TxTzRxRz", "withModules, TxRz", "withModules, TxTzRxRz", "2020, original c5, 20iter"]
names_c5 = ["2020 c5", "2019 c5", "2020, iter20 c5"] # number 4, 5 and 13
nums_c5 = [4, 5, 13]
# for 2019 data (up to #8)

fig1 = plt.figure()
fig1.suptitle("HighMomentum, chi2 / dofs vs. nTracks, 6500 events vs. 10000 events")
ax1 = fig1.add_subplot(121)

for run in range(len(names)):
  f = open('../output/alignment_runs/done_runs/HE_runs/{}/AlignmentResults/parsedlog.json'.format(run))
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

  line1 = ax1.plot(tracks, chi2dofs, label="{}".format(names[run]))
  ax1.legend()
  ax1.grid(True)
  ax1.set_xlabel("nTracks x 1000")
  ax1.set_ylabel("chi2 / dofs")
  f.close()
fig1.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_{}.pdf'.format(run))

# the same here for #9 to #12
ax2 = fig1.add_subplot(122)
for run in range(len(names_2020)):
  f = open('../output/alignment_runs/done_runs/HE_runs/{}/AlignmentResults/parsedlog.json'.format(run+9))
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
  print('run: ', run, ", min tracks: ", np.min(tracks), "max tracks: ", np.max(tracks), ", mean chi2: ", np.mean(chi2dofs))
  for i in range(len(tracks)):
    tracks[i] /= 1000

  line2 = ax2.plot(tracks, chi2dofs, label="{}".format(names_2020[run]))
  ax2.legend()
  ax2.grid(True)
  ax2.set_xlabel("nTracks x 1000")
  ax2.set_ylabel("chi2 / dofs")

  f.close()
fig1.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_{}.pdf'.format(run+9))
plt.show()
plt.clf()

# for the data up to #9

fig2 = plt.figure()
fig2.suptitle("HighMomentum, chi2 / dofs vs. iteration, 6500 events vs. 10000 events")
ax1 = fig2.add_subplot(121)

for run in range(len(names)):
  f = open('../output/alignment_runs/done_runs/HE_runs/{}/AlignmentResults/parsedlog.json'.format(run))
  data = json.load(f)

  vals, dofs = [], []
  for j in data['total_chi2_vals']:
    vals.append(j)
  for jj in data['total_chi2_dofs']:
    dofs.append(jj)

  vals, dofs = np.array(vals), np.array(dofs)
  chi2_per_dofs = vals / dofs

  # nTracks, also get these from json
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

  if run != 4:
    line1 = ax1.plot(np.linspace(1,10,10), chi2_per_dofs, label="{}".format(names[run]))
  if run == 4:
    line1 = ax1.plot(np.linspace(1,20,20), chi2_per_dofs, label="{}".format(names[run]))
  ax1.legend()
  ax1.grid(True)
  ax1.set_xlabel("iteration num")
  ax1.set_ylabel("chi2 / dofs")

  f.close()
fig2.savefig('outfiles/chi2_outfiles/chi2_vs_iternum_{}.pdf'.format(run))

# for numbers #9 to #12 for 2020 data
ax2 = fig2.add_subplot(1, 2, 2)
for run in range(len(names_2020)):
  f = open('../output/alignment_runs/done_runs/HE_runs/{}/AlignmentResults/parsedlog.json'.format(run+9))
  data = json.load(f)

  vals, dofs = [], []
  for j in data['total_chi2_vals']:
    vals.append(j)
  for jj in data['total_chi2_dofs']:
    dofs.append(jj)

  vals, dofs = np.array(vals), np.array(dofs)
  chi2_per_dofs = vals / dofs

  # nTracks, also get these from json
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
  
  if run != 4:
    line2 = ax2.plot(np.linspace(1,10,10), chi2_per_dofs, label="{}".format(names_2020[run]))
  if run == 4:
    line2 = ax2.plot(np.linspace(1,20,20), chi2_per_dofs, label="{}".format(names_2020[run]))
#  line2 = ax2.plot(np.linspace(1,10,10), chi2_per_dofs, label="{}".format(names_2020[run]))
  ax2.legend()
  ax2.grid(True)
  ax2.set_xlabel("iteration num")
  ax2.set_ylabel("chi2 / dofs")

  f.close()
fig2.savefig('outfiles/chi2_outfiles/chi2_vs_iternum_{}.pdf'.format(run))
plt.show()
plt.clf()

# for figure 3
fig3 = plt.figure()
fig3.suptitle("HighMomentum, chi2 / dofs vs. iteration, config5 for 2019 vs 2020")
ax1 = fig3.add_subplot(121)

for file_num in (nums_c5):
  f = open('../output/alignment_runs/done_runs/HE_runs/{}/AlignmentResults/parsedlog.json'.format(file_num))
  data = json.load(f)
  
  iternum = 0
  
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
  
  if file_num != 13:
    line1 = ax1.plot(np.linspace(1,10,10), chi2_per_dofs, label="{}".format(names_c5[iternum]))
  if file_num == 13:
    line1 = ax1.plot(np.linspace(1,20,20), chi2_per_dofs, label="{}".format(names_c5[iternum]))
  ax1.legend()
  ax1.grid(True)
  ax1.set_xlabel("iteration num")
  ax1.set_ylabel("chi2 / dofs")
  
  iternum += 1

  f.close()
fig3.savefig('outfiles/chi2_outfiles/chi2_vs_iternum_{}.pdf'.format(iternum))

ax2 = fig3.add_subplot(122)
for file_num in (nums_c5):
  f = open('../output/alignment_runs/done_runs/HE_runs/{}/AlignmentResults/parsedlog.json'.format(file_num))
  data = json.load(f)

  iternum = 0

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
  
  print('run: ', run, ", min tracks: ", np.min(tracks), "max tracks: ", np.max(tracks), ", chi2: ", chi2dofs)

  for i in range(len(tracks)):
    tracks[i] /= 1000

  line2 = ax2.plot(tracks, chi2dofs, label="{}".format(names_c5[iternum]))
  ax2.legend()
  ax2.grid(True)
  ax2.set_xlabel("nTracks x 1000")
  ax2.set_ylabel("chi2 / dofs")

  iternum += 1

  f.close()
fig3.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_{}.pdf'.format(iternum))
plt.show()
plt.clf()

######################################################

# number 13
name = "config 5, 2020, 20 iter"
f1 = open('../output/alignment_runs/done_runs/HE_runs/13/AlignmentResults/parsedlog.json')
data = json.load(f1)

vals, dofs = [], []
for j in data['total_chi2_vals']:
  vals.append(j)
for jj in data['total_chi2_dofs']:
  dofs.append(jj)

vals, dofs = np.array(vals), np.array(dofs)
chi2_per_dofs = vals / dofs

# nTracks, also get these from json
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

plt.plot(tracks, chi2dofs, label="chi2 / dofs")
plt.legend()
plt.grid()
plt.title("chi2 / dofs vs. nTracks: num13")
plt.xlabel("nTracks x 1000")
plt.ylabel("chi2 / dofs")
plt.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_13.pdf')
plt.show()
plt.clf()

plt.plot(np.linspace(1, 20, 20), nTracks, label="nTracks")
plt.plot(np.linspace(1, 20, 20), vals, label="total chi2 values")
plt.legend()
plt.grid()
plt.title("chi2 vs. nTracks: num13")
plt.xlabel("iteration")
plt.ylabel("tracks and chi2")
plt.savefig('outfiles/chi2_outfiles/tracks_chi2_vs_iter.pdf')
plt.show()
plt.clf()
