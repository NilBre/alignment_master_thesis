import re
import numpy as np
import matplotlib.pyplot as plt
import json

# loop over the 4 HE runs
names = ["onlyModules,2019,TxRz", "onlyModules,2020,TxRz"]
for run in range(len(names)):
  f = open('../output/alignment_runs/done_runs/HE_runs/{}/AlignmentResults/parsedlog.json'.format(run + 20))
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
  
  stations = ["FT/T1", "FT/T2", "FT/T3"]
  layers = ["LayerU", "LayerV", "LayerX1", "LayerX2"]  # different order
  quarters = ["Quarter(0|2)", "Quarter(1|3)"]
  T1T2_modules = ["Module0", "Module1", "Module2", "Module3", "Module4"]
  T3_modules = ["Module0", "Module1", "Module2", "Module3", "Module4", "Module5"]
  
  for s in stations:
    for l in layers:
      for q in quarters:
        if "T1" in s:
          modules = T1T2_modules
          for m in modules:
            section = s + l + q + m
            for i in data[section]['nTracks']:
              T1.append(float(i))
        if "T2" in s:
          modules = T1T2_modules
          for m in modules:
            section = s + l + q + m
            for ii in data[section]['nTracks']:
              T2.append(float(ii))
        if "T3" in s:
          modules = T3_modules
          for m in modules:
            section = s + l + q + m
            for iii in data[section]['nTracks']:
              T3.append(float(iii))

## reshape the lists to properly add each event num for each iteration
  new_T1 = np.reshape(T1, (40,10))
  new_T2 = np.reshape(T2, (40,10))
  new_T3 = np.reshape(T3, (48,10))
## sum over columns for overall tracks per iteration
  T1_iter_sum = np.sum(new_T1, axis=0)
  T2_iter_sum = np.sum(new_T2, axis=0)
  T3_iter_sum = np.sum(new_T3, axis=0)
  
  nTracks = np.array(T1_iter_sum) + np.array(T2_iter_sum) + np.array(T3_iter_sum)

  zipped_list = zip(nTracks, chi2_per_dofs)
  sorted_pairs = sorted(zipped_list)
  tuples = zip(*sorted_pairs)
  tracks, chi2dofs = [list(tuple) for tuple in tuples]

  for i in range(len(tracks)):
    tracks[i] /= 1000

  plt.plot(tracks, chi2dofs, label="{}".format(names[run]))
  plt.legend()
  plt.grid()
  plt.title("HighMomentum, chi2 / dofs vs. nTracks")
  plt.xlabel("nTracks x 1000")
  plt.ylabel("chi2 / dofs")

  f.close()

plt.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_{}.pdf'.format(run))
plt.show()
plt.clf()

################################################################

for run in range(len(names)):
  f = open('../output/alignment_runs/done_runs/HE_runs/{}/AlignmentResults/parsedlog.json'.format(run+20))
  data = json.load(f)

  vals, dofs = np.array(vals), np.array(dofs)
  chi2_per_dofs = vals / dofs

  # nTracks, also get these from json
  T1, T2, T3 = [], [], []

  stations = ["FT/T1", "FT/T2", "FT/T3"]
  layers = ["LayerU", "LayerV", "LayerX1", "LayerX2"]  # different order
  quarters = ["Quarter(0|2)", "Quarter(1|3)"]
  T1T2_modules = ["Module0", "Module1", "Module2", "Module3", "Module4"]
  T3_modules = ["Module0", "Module1", "Module2", "Module3", "Module4", "Module5"]

  for s in stations:
    for l in layers:
      for q in quarters:
        if "T1" in s:
          modules = T1T2_modules
          for m in modules:
            section = s + l + q + m
            for i in data[section]['nTracks']:
              T1.append(float(i))
        if "T2" in s:
          modules = T1T2_modules
          for m in modules:
            section = s + l + q + m
            for ii in data[section]['nTracks']:
              T2.append(float(ii))
        if "T3" in s:
          modules = T3_modules
          for m in modules:
            section = s + l + q + m
            for iii in data[section]['nTracks']:
              T3.append(float(iii))

## reshape the lists to properly add each event num for each iteration
  new_T1 = np.reshape(T1, (40,10))
  new_T2 = np.reshape(T2, (40,10))
  new_T3 = np.reshape(T3, (48,10))
## sum over columns for overall tracks per iteration
  T1_iter_sum = np.sum(new_T1, axis=0)
  T2_iter_sum = np.sum(new_T2, axis=0)
  T3_iter_sum = np.sum(new_T3, axis=0)

  nTracks = np.array(T1_iter_sum) + np.array(T2_iter_sum) + np.array(T3_iter_sum)

  plt.plot(np.linspace(1,10,10), nTracks, label="{}".format(names[run]))
  plt.legend()
  plt.grid()
  plt.title("iter vs. nTracks")
  plt.xlabel("iteration")
  plt.ylabel("nTracks x 1000")

  f.close()
plt.savefig('outfiles/chi2_outfiles/iter_vs_nTracks_{}.pdf'.format(run))
plt.show()
plt.clf()
