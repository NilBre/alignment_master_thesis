import re
import numpy as np
import matplotlib.pyplot as plt
import json
from collections import OrderedDict

#names = ["TxRz_2020_GL"]
# good long tracks from stack runs
names_GL = ["config5, TxRz", "config5, TxTzRxRz", "withModules, TxRz", "withModules,TxTzRxRz"]

# high momentum tracks
names_HE = ["config5, TxRz", "config5, TxTzRxRz", "withModules, TxRz", "withModules,TxTzRxRz"]

fig = plt.figure()
fig.suptitle("GoodLongTracks vs HighMomentumTTracks, 2020 data")
ax1 = fig.add_subplot(121)

# for goodlong tracks
for name in range(len(names_GL)):
  f = open("../../../../../../../../interactive_storage/nbreer/build_stack/stack/Alignment/output/no9to12/{}/AlignmentResults/parsedlog.json".format(name+9))
  data = json.load(f)

  vals, dofs = [], []
  for j in data['total_chi2_vals']:
    vals.append(j)
  for jj in data['total_chi2_dofs']:
    dofs.append(jj)
 
  vals, dofs = np.array(vals), np.array(dofs)
  chi2_per_dofs = vals / dofs

  ordered_chi2 = list(OrderedDict.fromkeys(chi2_per_dofs))
 
  T1, T2, T3 = [], [], []
  for i in data['FT/T1']['nTracks']:
    T1.append(float(i))
  for ii in data['FT/T2']['nTracks']:
    T2.append(float(ii))
  for iii in data['FT/T3']['nTracks']:
    T3.append(float(iii))

  nTracks = np.array(T1) + np.array(T2) + np.array(T3)

  zipped_list = zip(nTracks, ordered_chi2)
  sorted_pairs = sorted(zipped_list)
  tuples = zip(*sorted_pairs)
  tracks, chi2dofs = [list(tuple) for tuple in tuples]

  for i in range(len(tracks)):
    tracks[i] /= 1000

  ax1.plot(tracks, chi2dofs, label="{}".format(names_GL[name]))
  ax1.legend()
  ax1.grid(True)
  ax1.set_xlabel("nTracks x 1000")
  ax1.set_ylabel("chi2 / dofs")
  f.close()
fig.savefig('outfiles/chi2_vs_nTracks_{}.pdf'.format(name))

ax2 = fig.add_subplot(122)

# for HE tracks
for name in range(len(names_HE)):
  f = open("../../output/alignment_runs/done_runs/HE_runs/{}/AlignmentResults/parsedlog.json".format(name+9))
  data = json.load(f)

  vals, dofs = [], []
  for j in data['total_chi2_vals']:
    vals.append(j)
  for jj in data['total_chi2_dofs']:
    dofs.append(jj)
 
  vals, dofs = np.array(vals), np.array(dofs)
  chi2_per_dofs = vals / dofs

  ordered_chi2 = list(OrderedDict.fromkeys(chi2_per_dofs))
 
  T1, T2, T3 = [], [], []
  for i in data['FT/T1']['nTracks']:
    T1.append(float(i))
  for ii in data['FT/T2']['nTracks']:
    T2.append(float(ii))
  for iii in data['FT/T3']['nTracks']:
    T3.append(float(iii))

  nTracks = np.array(T1) + np.array(T2) + np.array(T3)

  zipped_list = zip(nTracks, ordered_chi2)
  sorted_pairs = sorted(zipped_list)
  tuples = zip(*sorted_pairs)
  tracks, chi2dofs = [list(tuple) for tuple in tuples]

  for i in range(len(tracks)):
    tracks[i] /= 1000

  ax2.plot(tracks, chi2dofs, label="{}".format(names_HE[name]))
  ax2.legend()
  ax2.grid(True)
  ax2.set_xlabel("nTracks x 1000")
  ax2.set_ylabel("chi2 / dofs")
  f.close()
fig.savefig('outfiles/chi2_vs_nTracks_{}.pdf'.format(name+9))
plt.show()
plt.clf()

###################### for iter vs chi2 #########################

fig2 = plt.figure()
fig2.suptitle("GoodLongTracks vs HighMomentumTTracks, 2020 data")
ax1 = fig2.add_subplot(121)

for name in range(len(names_GL)):
  f = open("../../../../../../../../interactive_storage/nbreer/build_stack/stack/Alignment/output/no9to12/{}/AlignmentResults/parsedlog.json".format(name+9))
  data = json.load(f)

  vals, dofs = [], []
  for j in data['total_chi2_vals']:
    vals.append(j)
  for jj in data['total_chi2_dofs']:
    dofs.append(jj)

  vals, dofs = np.array(vals), np.array(dofs)
  chi2_per_dofs = vals / dofs

  ordered_chi2 = list(OrderedDict.fromkeys(chi2_per_dofs))

  ax1.plot(np.linspace(1,10,10), ordered_chi2, label="{}".format(names_GL[name]))
  ax1.legend()
  ax1.grid(True)
  ax1.set_xlabel("iteration num")
  ax1.set_ylabel("chi2 / dofs")

  f.close()
fig2.savefig('outfiles/chi2_vs_iternum_{}.pdf'.format(name+9))

ax2 = fig2.add_subplot(122)

for name in range(len(names_HE)):
  f = open("../../output/alignment_runs/done_runs/HE_runs/{}/AlignmentResults/parsedlog.json".format(name+9))
  data = json.load(f)

  vals, dofs = [], []
  for j in data['total_chi2_vals']:
    vals.append(j)
  for jj in data['total_chi2_dofs']:
    dofs.append(jj)

  vals, dofs = np.array(vals), np.array(dofs)
  chi2_per_dofs = vals / dofs

  ordered_chi2 = list(OrderedDict.fromkeys(chi2_per_dofs))

  ax2.plot(np.linspace(1,10,10), ordered_chi2, label="{}".format(names_GL[name]))
  ax2.legend()
  ax2.grid(True)
  ax2.set_xlabel("iteration num")
  ax2.set_ylabel("chi2 / dofs")

  f.close()
fig2.savefig('outfiles/chi2_vs_iternum_{}.pdf'.format(name+9))

plt.show()
plt.clf()

############################### names GL in seperate plots ########################
names = ["config5, TxRz", "config5, TxTzRxRz", "withModules, TxRz", "withModules,TxTzRxRz"]

fig3 = plt.figure()
fig3.suptitle("GoodLongTracks 2020 data")

for name in range(len(names)):
  f = open("../../../../../../../../interactive_storage/nbreer/build_stack/stack/Alignment/output/no9to12/{}/AlignmentResults/parsedlog.json".format(name+9))
  data = json.load(f)

  vals, dofs = [], []
  for j in data['total_chi2_vals']:
    vals.append(j)
  for jj in data['total_chi2_dofs']:
    dofs.append(jj)

  vals, dofs = np.array(vals), np.array(dofs)
  chi2_per_dofs = vals / dofs

  ordered_chi2 = list(OrderedDict.fromkeys(chi2_per_dofs))

  T1, T2, T3 = [], [], []
  for i in data['FT/T1']['nTracks']:
    T1.append(float(i))
  for ii in data['FT/T2']['nTracks']:
    T2.append(float(ii))
  for iii in data['FT/T3']['nTracks']:
    T3.append(float(iii))

  nTracks = np.array(T1) + np.array(T2) + np.array(T3)

  zipped_list = zip(nTracks, ordered_chi2)
  sorted_pairs = sorted(zipped_list)
  tuples = zip(*sorted_pairs)
  tracks, chi2dofs = [list(tuple) for tuple in tuples]

  for i in range(len(tracks)):
    tracks[i] /= 1000

  if name == 0:
    ax1 = fig3.add_subplot(221)
    ax1.plot(tracks, chi2dofs, label="{}".format(names_GL[name]))
    ax1.legend()
    ax1.grid(True)
    ax1.set_xlabel("nTracks x 1000")
    ax1.set_ylabel("chi2 / dofs")
    f.close()
  if name == 1:
    ax2 = fig3.add_subplot(222)
    ax2.plot(tracks, chi2dofs, label="{}".format(names_GL[name]))
    ax2.legend()
    ax2.grid(True)
    ax2.set_xlabel("nTracks x 1000")
    ax2.set_ylabel("chi2 / dofs")
    f.close()
  if name == 2:
    ax3 = fig3.add_subplot(223)
    ax3.plot(tracks, chi2dofs, label="{}".format(names_GL[name]))
    ax3.legend()
    ax3.grid(True)
    ax3.set_xlabel("nTracks x 1000")
    ax3.set_ylabel("chi2 / dofs")
    f.close()
  if name == 3:
    ax4 = fig3.add_subplot(224)
    ax4.plot(tracks, chi2dofs, label="{}".format(names_GL[name]))
    ax4.legend()
    ax4.grid(True)
    ax4.set_xlabel("nTracks x 1000")
    ax4.set_ylabel("chi2 / dofs")
    f.close()
fig.savefig('outfiles/chi2_vs_nTracks_{}.pdf'.format(name))
plt.show()
plt.clf()
