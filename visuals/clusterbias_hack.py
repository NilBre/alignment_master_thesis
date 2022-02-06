
import re
import numpy as np
import matplotlib.pyplot as plt
import json

##############################################################
names_stack = ["modules_low_lumi", "modules_normal_lumi"]
names_cluster_bias = ["modules_low", "modules_normal"]
names_constr = ["without_1st_constr", "without_2nd_constr"]

# define data for both low and normal lumi
f1 = open('../output/stack_alignment/modules_low_lumi/AlignmentResults/parsedlog.json')
data1 = json.load(f1)

f2 = open('../output/stack_alignment/modules_normal_lumi/AlignmentResults/parsedlog.json')
data2 = json.load(f2)

f3 = open('../output/stack_alignment/clusterbias_hack/modules_low/AlignmentResults/parsedlog.json')
data3 = json.load(f3)

f4 = open('../output/stack_alignment/clusterbias_hack/modules_normal/AlignmentResults/parsedlog.json')
data4 = json.load(f4)

f5 = open('../output/stack_alignment/clusterbias_hack/without_1st_constr/AlignmentResults/parsedlog.json')
data5 = json.load(f5)

f6 = open('../output/stack_alignment/clusterbias_hack/without_2nd_constr/AlignmentResults/parsedlog.json')
data6 = json.load(f6)

########## clusterbias with misalignment + low lumi #########
f7 = open("../output/misalignment_runs/clusterbias/0/AlignmentResults/parsedlog.json")
data7 = json.load(f7)

f8 = open("../output/misalignment_runs/clusterbias/1/AlignmentResults/parsedlog.json")
data8 = json.load(f8)

f9 = open("../output/misalignment_runs/clusterbias/2/AlignmentResults/parsedlog.json")
data9 = json.load(f9)

f10 = open("../output/misalignment_runs/clusterbias/3/AlignmentResults/parsedlog.json")
data10 = json.load(f10)

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

vals_cbl, dofs_cbl = [], []
for j in data3['total_chi2_vals']:
  vals_cbl.append(j)
for jj in data3['total_chi2_dofs']:
  dofs_cbl.append(jj)

vals_cbn, dofs_cbn = [], []
for j in data4['total_chi2_vals']:
  vals_cbn.append(j)
for jj in data4['total_chi2_dofs']:
  dofs_cbn.append(jj)

vals_1, dofs_1 = [], []
for j in data5['total_chi2_vals']:
  vals_1.append(j)
for jj in data5['total_chi2_dofs']:
  dofs_1.append(jj)

vals_2, dofs_2 = [], []
for j in data6['total_chi2_vals']:
  vals_2.append(j)
for jj in data6['total_chi2_dofs']:
  dofs_2.append(jj)

######## misalign
vals_mis0, dofs_mis0 = [], []
for i in data7['total_chi2_vals']:
  vals_mis0.append(i)
for ii in data7['total_chi2_dofs']:
  dofs_mis0.append(ii)

vals_mis1, dofs_mis1 = [], []
for i in data8['total_chi2_vals']:
  vals_mis1.append(i)
for ii in data8['total_chi2_dofs']:
  dofs_mis1.append(ii)

vals_mis2, dofs_mis2 = [], []
for i in data9['total_chi2_vals']:
  vals_mis2.append(i)
for ii in data9['total_chi2_dofs']:
  dofs_mis2.append(ii)

vals_mis3, dofs_mis3 = [], []
for i in data10['total_chi2_vals']:
  vals_mis3.append(i)
for ii in data10['total_chi2_dofs']:
  dofs_mis3.append(ii)

########

vals_l, dofs_l = np.array(vals_l), np.array(dofs_l)
chi2_per_dofs_l = vals_l / dofs_l

vals_n, dofs_n = np.array(vals_n), np.array(dofs_n)
chi2_per_dofs_n = vals_n / dofs_n

vals_cbl, dofs_cbl = np.array(vals_cbl), np.array(dofs_cbl)
chi2_per_dofs_cbl = vals_cbl / dofs_cbl

vals_cbn, dofs_cbn = np.array(vals_cbn), np.array(dofs_cbn)
chi2_per_dofs_cbn = vals_cbn / dofs_cbn

vals_1, dofs_1 = np.array(vals_1), np.array(dofs_1)
chi2_per_dofs_1 = vals_1 / dofs_1

vals_2, dofs_2 = np.array(vals_2), np.array(dofs_2)
chi2_per_dofs_2 = vals_2 / dofs_2

### misalign ###
vals_mis0, dofs_mis0 = np.array(vals_mis0), np.array(dofs_mis0)
chi2_per_dofs_mis0 = vals_mis0 / dofs_mis0

vals_mis1, dofs_mis1 = np.array(vals_mis1), np.array(dofs_mis1)
chi2_per_dofs_mis1 = vals_mis1 / dofs_mis1

vals_mis0, dofs_mis2 = np.array(vals_mis2), np.array(dofs_mis2)
chi2_per_dofs_mis2 = vals_mis2 / dofs_mis2

vals_mis0, dofs_mis3 = np.array(vals_mis3), np.array(dofs_mis3)
chi2_per_dofs_mis3 = vals_mis3 / dofs_mis3

# nTracks, also get these from json
T1_l, T2_l, T3_l = [], [], []
T1_n, T2_n, T3_n = [], [], []

T1_cbl, T2_cbl, T3_cbl = [], [], []
T1_cbn, T2_cbn, T3_cbn = [], [], []

T1_1, T2_1, T3_1 = [], [], []
T1_2, T2_2, T3_2 = [], [], []

T1_mis0, T2_mis0, T3_mis0 = [], [], []
T1_mis1, T2_mis1, T3_mis1 = [], [], []
T1_mis2, T2_mis2, T3_mis2 = [], [], []
T1_mis3, T2_mis3, T3_mis3 = [], [], []

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
          for iii in data3[section]['nTracks']:
            T1_cbl.append(float(iii))
          for iv in data4[section]['nTracks']:
            T1_cbn.append(float(iv))
          for m in data7[section]['nTracks']:
            T1_mis0.append(float(m))
          for mm in data8[section]['nTracks']:
            T1_mis1.append(float(mm))
          for mmm in data9[section]['nTracks']:
            T1_mis2.append(float(mmm))
          for im in data10[section]['nTracks']:
            T1_mis3.append(float(im))
      if "T2" in s:
        modules = T1T2_modules
        for m in modules:
          section = s + l + q + m
          for j in data1[section]['nTracks']:
            T2_l.append(float(j))
          for jj in data2[section]['nTracks']:
            T2_n.append(float(jj))
          for jjj in data3[section]['nTracks']:
            T2_cbl.append(float(jjj))
          for jv in data4[section]['nTracks']:
            T2_cbn.append(float(jv))
          for o in data7[section]['nTracks']:
            T2_mis0.append(float(o))
          for oo in data8[section]['nTracks']:
            T2_mis1.append(float(oo))
          for ooo in data9[section]['nTracks']:
            T2_mis2.append(float(ooo))
          for io in data10[section]['nTracks']:
            T2_mis3.append(float(io))
      if "T3" in s:
        modules = T3_modules
        for m in modules:
          section = s + l + q + m
          for k in data1[section]['nTracks']:
            T3_l.append(float(k))
          for kk in data2[section]['nTracks']:
            T3_n.append(float(kk))
          for kkk in data3[section]['nTracks']:
            T3_cbl.append(float(kkk))
          for kv in data4[section]['nTracks']:
            T3_cbn.append(float(kv))
          for p in data7[section]['nTracks']:
            T3_mis0.append(float(p))
          for pp in data8[section]['nTracks']:
            T3_mis1.append(float(pp))
          for ppp in data9[section]['nTracks']:
            T3_mis2.append(float(ppp))
          for ip in data10[section]['nTracks']:
            T3_mis3.append(float(ip))

## reshape the lists to properly add each event num for each iteration
new_T1_l = np.reshape(T1_l, (40,10))
new_T2_l = np.reshape(T2_l, (40,10))
new_T3_l = np.reshape(T3_l, (48,10))

new_T1_n = np.reshape(T1_n, (40,10))
new_T2_n = np.reshape(T2_n, (40,10))
new_T3_n = np.reshape(T3_n, (48,10))

new_T1_cbl = np.reshape(T1_cbl, (40,10))
new_T2_cbl = np.reshape(T2_cbl, (40,10))
new_T3_cbl = np.reshape(T3_cbl, (48,10))

new_T1_cbn = np.reshape(T1_cbn, (40,10))
new_T2_cbn = np.reshape(T2_cbn, (40,10))
new_T3_cbn = np.reshape(T3_cbn, (48,10))

######## misalign #########
new_T1_mis0 = np.reshape(T1_mis0, (40,10))
new_T2_mis0 = np.reshape(T2_mis0, (40,10))
new_T3_mis0 = np.reshape(T3_mis0, (48,10))

new_T1_mis1 = np.reshape(T1_mis1, (40,10))
new_T2_mis1 = np.reshape(T2_mis1, (40,10))
new_T3_mis1 = np.reshape(T3_mis1, (48,10))

new_T1_mis2 = np.reshape(T1_mis2, (40,10))
new_T2_mis2 = np.reshape(T2_mis2, (40,10))
new_T3_mis2 = np.reshape(T3_mis2, (48,10))

new_T1_mis3 = np.reshape(T1_mis3, (40,10))
new_T2_mis3 = np.reshape(T2_mis3, (40,10))
new_T3_mis3 = np.reshape(T3_mis3, (48,10))
## sum over columns for overall tracks per iteration
T1_iter_sum_l = np.sum(new_T1_l, axis=0)
T2_iter_sum_l = np.sum(new_T2_l, axis=0)
T3_iter_sum_l = np.sum(new_T3_l, axis=0)

T1_iter_sum_n = np.sum(new_T1_n, axis=0)
T2_iter_sum_n = np.sum(new_T2_n, axis=0)
T3_iter_sum_n = np.sum(new_T3_n, axis=0)

T1_iter_sum_cbl = np.sum(new_T1_cbl, axis=0)
T2_iter_sum_cbl = np.sum(new_T2_cbl, axis=0)
T3_iter_sum_cbl = np.sum(new_T3_cbl, axis=0)

T1_iter_sum_cbn = np.sum(new_T1_cbn, axis=0)
T2_iter_sum_cbn = np.sum(new_T2_cbn, axis=0)
T3_iter_sum_cbn = np.sum(new_T3_cbn, axis=0)

### with misalignment
T1_iter_sum_mis0 = np.sum(new_T1_mis0, axis=0)
T2_iter_sum_mis0 = np.sum(new_T2_mis0, axis=0)
T3_iter_sum_mis0 = np.sum(new_T3_mis0, axis=0)

T1_iter_sum_mis1 = np.sum(new_T1_mis1, axis=0)
T2_iter_sum_mis1 = np.sum(new_T2_mis1, axis=0)
T3_iter_sum_mis1 = np.sum(new_T3_mis1, axis=0)

T1_iter_sum_mis2 = np.sum(new_T1_mis2, axis=0)
T2_iter_sum_mis2 = np.sum(new_T2_mis2, axis=0)
T3_iter_sum_mis2 = np.sum(new_T3_mis2, axis=0)

T1_iter_sum_mis3 = np.sum(new_T1_mis3, axis=0)
T2_iter_sum_mis3 = np.sum(new_T2_mis3, axis=0)
T3_iter_sum_mis3 = np.sum(new_T3_mis3, axis=0)
####################################################
# when i use stations + layers
for i in data5['FT/T1']['nTracks']:
  T1_1.append(float(i))
for ii in data5['FT/T2']['nTracks']:
  T2_1.append(float(ii))
for iii in data5['FT/T3']['nTracks']:
  T3_1.append(float(iii))

for j in data6['FT/T1']['nTracks']:
  T1_2.append(float(j))
for jj in data6['FT/T2']['nTracks']:
  T2_2.append(float(jj))
for jjj in data6['FT/T3']['nTracks']:
  T3_2.append(float(jjj))

####################################################
nTracks_l = np.array(T1_iter_sum_l) + np.array(T2_iter_sum_l) + np.array(T3_iter_sum_l)

nTracks_n = np.array(T1_iter_sum_n) + np.array(T2_iter_sum_n) + np.array(T3_iter_sum_n)

nTracks_cbl = np.array(T1_iter_sum_cbl) + np.array(T2_iter_sum_cbl) + np.array(T3_iter_sum_cbl)

nTracks_cbn = np.array(T1_iter_sum_cbn) + np.array(T2_iter_sum_cbn) + np.array(T3_iter_sum_cbn)

nTracks_1 = np.array(T1_1) + np.array(T2_1) + np.array(T3_1)

nTracks_2 = np.array(T1_2) + np.array(T2_2) + np.array(T3_2)

# misaligned
nTracks_mis0 = np.array(T1_iter_sum_mis0) + np.array(T2_iter_sum_mis0) + np.array(T3_iter_sum_mis0)
nTracks_mis1 = np.array(T1_iter_sum_mis1) + np.array(T2_iter_sum_mis1) + np.array(T3_iter_sum_mis1)
nTracks_mis2 = np.array(T1_iter_sum_mis2) + np.array(T2_iter_sum_mis2) + np.array(T3_iter_sum_mis2)
nTracks_mis3 = np.array(T1_iter_sum_mis3) + np.array(T2_iter_sum_mis3) + np.array(T3_iter_sum_mis3)
####################################################
# removing constraints sample
zipped_list_1 = zip(nTracks_1, chi2_per_dofs_1)
sorted_pairs_1 = sorted(zipped_list_1)
tuples_1 = zip(*sorted_pairs_1)
tracks_1, chi2dofs_1= [list(tuple_1) for tuple_1 in tuples_1]

for i in range(len(tracks_1)):
  tracks_1[i] /= 1000

zipped_list_2 = zip(nTracks_2, chi2_per_dofs_2)
sorted_pairs_2 = sorted(zipped_list_2)
tuples_2 = zip(*sorted_pairs_2)
tracks_2, chi2dofs_2= [list(tuple_2) for tuple_2 in tuples_2]

for i in range(len(tracks_2)):
  tracks_2[i] /= 1000

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

###
zipped_list_cbl = zip(nTracks_cbl, chi2_per_dofs_cbl)
sorted_pairs_cbl = sorted(zipped_list_cbl)
tuples_cbl = zip(*sorted_pairs_cbl)
tracks_cbl, chi2dofs_cbl= [list(tuple_cbl) for tuple_cbl in tuples_cbl]

for i in range(len(tracks_cbl)):
  tracks_cbl[i] /= 1000

###
zipped_list_cbn = zip(nTracks_n, chi2_per_dofs_cbn)
sorted_pairs_cbn = sorted(zipped_list_cbn)
tuples_cbn = zip(*sorted_pairs_cbn)
tracks_cbn, chi2dofs_cbn = [list(tuple_cbn) for tuple_cbn in tuples_cbn]

for i in range(len(tracks_cbn)):
  tracks_cbn[i] /= 1000

### misaligned ###
zipped_list_mis0 = zip(nTracks_mis0, chi2_per_dofs_mis0)
sorted_pairs_mis0 = sorted(zipped_list_mis0)
tuples_mis0 = zip(*sorted_pairs_mis0)
tracks_mis0, chi2dofs_mis0 = [list(tuple_mis0) for tuple_mis0 in tuples_mis0]

for i in range(len(tracks_mis0)):
  tracks_mis0[i] /= 1000

###
zipped_list_mis1 = zip(nTracks_mis1, chi2_per_dofs_mis1)
sorted_pairs_mis1 = sorted(zipped_list_mis1)
tuples_mis1 = zip(*sorted_pairs_mis1)
tracks_mis1, chi2dofs_mis1 = [list(tuple_mis1) for tuple_mis1 in tuples_mis1]

for i in range(len(tracks_mis1)):
  tracks_mis1[i] /= 1000

###
zipped_list_mis2 = zip(nTracks_mis2, chi2_per_dofs_mis2)
sorted_pairs_mis2 = sorted(zipped_list_mis2)
tuples_mis2 = zip(*sorted_pairs_mis2)
tracks_mis2, chi2dofs_mis2 = [list(tuple_mis2) for tuple_mis2 in tuples_mis2]

for i in range(len(tracks_mis2)):
  tracks_mis2[i] /= 1000

###
zipped_list_mis3 = zip(nTracks_mis3, chi2_per_dofs_mis3)
sorted_pairs_mis3 = sorted(zipped_list_mis3)
tuples_mis3 = zip(*sorted_pairs_mis3)
tracks_mis3, chi2dofs_mis3 = [list(tuple_mis3) for tuple_mis3 in tuples_mis3]

for i in range(len(tracks_mis3)):
  tracks_mis3[i] /= 1000

#######################
fig,ax = plt.subplots()
ax.plot(tracks_l, chi2dofs_l, color="red", label="low lumi")
ax.legend()
ax.grid(True)
ax.set_title("GoodLongTracks, ModulesOnly, config5 constraints")
ax.set_xlabel("nTracks x 1000")
ax.set_ylabel("chi2 / dofs", color="red")

ax2 = ax.twinx()
ax2.plot(tracks_l, chi2dofs_n, color="blue", label="normal lumi")
ax2.legend()
ax2.grid(True)
ax2.set_ylabel("chi2 / dofs", color="blue")

f1.close()
f2.close()

fig.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_low.pdf', bbox_inches="tight")
plt.show()
plt.clf()

########

fig3,ax3 = plt.subplots()

reverse_cbl = np.sort(chi2dofs_cbl)[::-1]
reverse_cbn = np.sort(chi2dofs_cbn)[::-1]

ax3.plot(tracks_cbl, reverse_cbl, color="red", label="low lumi")
ax3.legend()
ax3.grid(True)
ax3.set_title("GoodLongTracks, ModulesOnly, config5 constraints, with clusterbias hack")
ax3.set_xlabel("nTracks x 1000")
ax3.set_ylabel("chi2 / dofs", color="red")

ax4 = ax3.twinx()
ax4.plot(tracks_cbl, reverse_cbn, color="blue", label="normal lumi")
ax4.legend()
ax4.grid(True)
ax4.set_ylabel("chi2 / dofs", color="blue")

f3.close()
f4.close()

fig3.savefig('outfiles/chi2_outfiles/chi2_vs_nTracks_low.pdf', bbox_inches="tight")
plt.show()
plt.clf()

#######################
# TODO:
# 1. plot the samples without first and second constraint
# 2. plot for each of the misaligned samples
# 3. one plot compare all misaligned samples
# 4. plot misaligned low vs normal lumi for every misalignment
#######################

# fix more scaling issues

#reverse_1 = np.sort(chi2dofs_1)[::-1]
#reverse_2 = np.sort(chi2dofs_2)[::-1]

fig5,ax5 = plt.subplots()
ax5.plot(tracks_1, chi2dofs_1, color="red", label="without constraint 1")
ax5.legend()
ax5.grid(True)
ax5.set_title("normal clusterbias vs without 1st constraint, low lumi")
ax5.set_xlabel("nTracks x 1000")
ax5.set_ylabel("chi2 / dofs", color="red")

ax6 = ax5.twinx()
ax6.plot(tracks_1, chi2dofs_2, color="blue", label="without constraint 2")
ax6.legend()
ax6.grid(True)
ax6.set_ylabel("chi2 / dofs", color="blue")

f5.close()
f6.close()

fig5.savefig('outfiles/chi2_outfiles/chi2dofs_removed_constraint1_low.pdf', bbox_inches="tight")
plt.show()
plt.clf()

#####
print("chi2dofs_100mu", chi2dofs_mis2)
fig7,ax7 = plt.subplots()
ax7.plot(tracks_mis2, chi2dofs_mis2, color="red", label="misalign 1mu")
ax7.legend()
ax7.grid(True)
ax7.set_title("CB, translation misalignment, low lumi")
ax7.set_xlabel("nTracks x 1000")
ax7.set_ylabel("chi2 / dofs", color="red")

ax8 = ax7.twinx()
ax8.plot(tracks_mis2, chi2dofs_mis3, color="blue", label="misalign 0.1mu")
ax8.legend()
ax8.grid(True)
ax8.set_ylabel("chi2 / dofs", color="blue")

f7.close()
f8.close()

fig7.savefig('outfiles/chi2_outfiles/chi2dofs_removed_constraint2_low.pdf', bbox_inches="tight")
plt.show()
plt.clf()

###
