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

##############################################################
##############################################################
#names_stack = ["modules_low_lumi", "modules_normal_lumi"]

# with clusterbias hack
names_stack = ["modules_low", "modules_normal"]

# define data for both low and normal lumi
#f1 = open('../output/stack_alignment/modules_low_lumi/AlignmentResults/parsedlog.json')
f1 = open('../output/stack_alignment/clusterbias_hack/modules_low/AlignmentResults/parsedlog.json')
data1 = json.load(f1)

#f2 = open('../output/stack_alignment/modules_normal_lumi/AlignmentResults/parsedlog.json')
f2 = open('../output/stack_alignment/clusterbias_hack/modules_normal/AlignmentResults/parsedlog.json')
data2 = json.load(f2)

fig,ax = plt.subplots()

#store dofs and vals for both low and normal
vals_l, dofs_l = [], []
for j in data1['total_chi2_vals']:
  vals_l.append(j)
for jj in data1['total_chi2_dofs']:
  dofs_l.append(jj)

vals_n, dofs_n = [], []
for j in data2['total_chi2_vals']:
  vals_n.append(j)
for jj in data2['total_chi2_dofs']:
  dofs_n.append(jj)

vals_l, dofs_l = np.array(vals_l), np.array(dofs_l)
chi2_per_dofs_l = vals_l / dofs_l

vals_n, dofs_n = np.array(vals_n), np.array(dofs_n)
chi2_per_dofs_n = vals_n / dofs_n

# nTracks, also get these from json
T1_l, T2_l, T3_l = [], [], []
T1_n, T2_n, T3_n = [], [], []

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
          for i in data1[section]['nTracks']:
            T1_l.append(float(i))
          for ii in data2[section]['nTracks']:
            T1_n.append(float(ii))
      if "T2" in s:
        modules = T1T2_modules
        for m in modules:
          section = s + l + q + m
          for j in data1[section]['nTracks']:
            T2_l.append(float(j))
          for jj in data2[section]['nTracks']:
            T2_n.append(float(jj)) 
      if "T3" in s:
        modules = T3_modules
        for m in modules:
          section = s + l + q + m
          for k in data1[section]['nTracks']:
            T3_l.append(float(k))
          for kk in data2[section]['nTracks']:
            T3_n.append(float(kk))

## reshape the lists to properly add each event num for each iteration
new_T1_l = np.reshape(T1_l, (40,10))
new_T2_l = np.reshape(T2_l, (40,10))
new_T3_l = np.reshape(T3_l, (48,10))

new_T1_n = np.reshape(T1_n, (40,10))
new_T2_n = np.reshape(T2_n, (40,10))
new_T3_n = np.reshape(T3_n, (48,10))
## sum over columns for overall tracks per iteration
T1_iter_sum_l = np.sum(new_T1_l, axis=0)
T2_iter_sum_l = np.sum(new_T2_l, axis=0)
T3_iter_sum_l = np.sum(new_T3_l, axis=0)

T1_iter_sum_n = np.sum(new_T1_n, axis=0)
T2_iter_sum_n = np.sum(new_T2_n, axis=0)
T3_iter_sum_n = np.sum(new_T3_n, axis=0)
####################################################
nTracks_l = np.array(T1_iter_sum_l) + np.array(T2_iter_sum_l) + np.array(T3_iter_sum_l)

nTracks_n = np.array(T1_iter_sum_n) + np.array(T2_iter_sum_n) + np.array(T3_iter_sum_n)
####################################################
# sort data for low and normal
zipped_list_l = zip(nTracks_l, chi2_per_dofs_l)
sorted_pairs_l = sorted(zipped_list_l)
tuples_l = zip(*sorted_pairs_l)
tracks_l, chi2dofs_l= [list(tuple_l) for tuple_l in tuples_l]

for i in range(len(tracks_l)):
  tracks_l[i] /= 1000

###
zipped_list_n = zip(nTracks_n, chi2_per_dofs_n)
sorted_pairs_n = sorted(zipped_list_n)
tuples_n = zip(*sorted_pairs_n)
tracks_n, chi2dofs_n = [list(tuple_n) for tuple_n in tuples_n]

for i in range(len(tracks_n)):
  tracks_n[i] /= 1000

#######################

chi2dofs_l = [1.3062874643810793, 
1.3062595317446717, 
1.3062796533186278, 
1.3062797601642093, 
1.3062693370443663, 
1.3062708328698354, 
1.3062662385596497, 
1.3062708328698354, 
1.3062610032223225, 
1.3062736489142226]

chi2dofs_n = [1.2855496836208735, 
1.2855524096392938, 
1.285551405315373, 
1.2855501140440662, 
1.2855345361640502, 
1.285552092214044, 
1.2855505444675468, 
1.285552092214044, 
1.2856138702993174, 
1.2855608136144963]

chi2dofs_l_rd = [round(i,6) for i in chi2dofs_l]
chi2dofs_n_rd = [round(i,6) for i in chi2dofs_n]

print("chi2_l mean: ", np.mean(chi2dofs_l))
print("chi2_n mean: ", np.mean(chi2dofs_n))

print("delta min max low: ", abs(np.max(chi2dofs_l_rd) - np.min(chi2dofs_l_rd)))
print("delta min max normal: ", abs(np.max(chi2dofs_n_rd) - np.min(chi2dofs_n_rd)))

ax.plot(tracks_l, chi2dofs_l_rd, color="red", label="low lumi")
ax.legend()
ax.grid(True)
ax.set_title("GoodLongTracks, ModulesOnly, config5 constraints")
ax.set_xlabel("nTracks x 1000")
ax.set_ylabel("chi2 / dofs", color="red")
ax.set_ylim(1.306255, 1.306288)

ax2 = ax.twinx()
ax2.plot(tracks_l, chi2dofs_n_rd, color="blue", label="normal lumi")
ax2.legend()
ax2.grid(True)
ax2.set_ylabel("chi2 / dofs", color="blue")
ax2.set_ylim(1.285533, 1.285614)

f1.close()
f2.close()

fig.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_low.pdf', bbox_inches="tight")
plt.show()
plt.clf()

######
#f3 = open('../output/stack_alignment/modules_low_lumi/AlignmentResults/parsedlog.json')
f3 = open('../output/stack_alignment/clusterbias_hack/modules_low/AlignmentResults/parsedlog.json')
data3 = json.load(f3)

#f4 = open('../output/stack_alignment/modules_normal_lumi/AlignmentResults/parsedlog.json')
f4 = open('../output/stack_alignment/clusterbias_hack/modules_normal/AlignmentResults/parsedlog.json')
data4 = json.load(f4)

fig1,ax3 = plt.subplots()

#store dofs and vals for both low and normal
vals_l, dofs_l = [], []
for j in data3['total_chi2_vals']:
  vals_l.append(j)
for jj in data3['total_chi2_dofs']:
  dofs_l.append(jj)

vals_n, dofs_n = [], []
for j in data2['total_chi2_vals']:
  vals_n.append(j)
for jj in data2['total_chi2_dofs']:
  dofs_n.append(jj)

vals_l, dofs_l = np.array(vals_l), np.array(dofs_l)
chi2_per_dofs_l = vals_l / dofs_l

vals_n, dofs_n = np.array(vals_n), np.array(dofs_n)
chi2_per_dofs_n = vals_n / dofs_n

# nTracks, also get these from json
T1_l, T2_l, T3_l = [], [], []
T1_n, T2_n, T3_n = [], [], []

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
          for i in data3[section]['nTracks']:
            T1_l.append(float(i))
          for ii in data4[section]['nTracks']:
            T1_n.append(float(ii))
      if "T2" in s:
        modules = T1T2_modules
        for m in modules:
          section = s + l + q + m
          for j in data3[section]['nTracks']:
            T2_l.append(float(j))
          for jj in data4[section]['nTracks']:
            T2_n.append(float(jj))
      if "T3" in s:
        modules = T3_modules
        for m in modules:
          section = s + l + q + m
          for k in data3[section]['nTracks']:
            T3_l.append(float(k))
          for kk in data4[section]['nTracks']:
            T3_n.append(float(kk))

## reshape the lists to properly add each event num for each iteration
new_T1_l = np.reshape(T1_l, (40,10))
new_T2_l = np.reshape(T2_l, (40,10))
new_T3_l = np.reshape(T3_l, (48,10))

new_T1_n = np.reshape(T1_n, (40,10))
new_T2_n = np.reshape(T2_n, (40,10))
new_T3_n = np.reshape(T3_n, (48,10))
## sum over columns for overall tracks per iteration
T1_iter_sum_l = np.sum(new_T1_l, axis=0)
T2_iter_sum_l = np.sum(new_T2_l, axis=0)
T3_iter_sum_l = np.sum(new_T3_l, axis=0)

T1_iter_sum_n = np.sum(new_T1_n, axis=0)
T2_iter_sum_n = np.sum(new_T2_n, axis=0)
T3_iter_sum_n = np.sum(new_T3_n, axis=0)
####################################################
nTracks_l = np.array(T1_iter_sum_l) + np.array(T2_iter_sum_l) + np.array(T3_iter_sum_l)

nTracks_n = np.array(T1_iter_sum_n) + np.array(T2_iter_sum_n) + np.array(T3_iter_sum_n)
####################################################
ax3.plot(np.linspace(0,9,10), nTracks_l, color="red", label="low lumi")
ax3.legend()
ax3.grid(True)
ax3.set_title("GoodLongTracks, ModulesOnly, config5 constraints")
ax3.set_xlabel("iterations")
ax3.set_ylabel("nTracks, nu=3.8", color="red")

ax4 = ax3.twinx()
ax4.plot(np.linspace(0,9,10), nTracks_n, color="blue", label="normal lumi")
ax4.legend()
ax4.grid(True)
ax4.set_ylabel("nTracks, nu=7.6", color="blue")

f3.close()
f4.close()

fig1.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_low.pdf', bbox_inches="tight")
plt.show()
plt.clf()
