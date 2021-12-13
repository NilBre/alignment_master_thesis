import matplotlib.pyplot as plt

#for recolor
import matplotlib.colors as mplcolor
import colorsys as colorsys

import json
import numpy as np

labels = ["Tx","Ty","Tz"]
positions = ["x_global","y_global","z_global"]

expectations = [0,0,0,0,0,0]

plotOn = True

title_t = "Layer translation averaged"
title_r = "Layer rotation averaged"

def layer_to_int(inputstring):
    if "LayerX1" in inputstring:
        return 0
    elif "LayerU" in inputstring:
        return 1
    elif "LayerV" in inputstring:
        return 2
    elif "LayerX2" in inputstring:
        return 3
    else:
        return -1

def halfstation_to_int(inputstring):
    if "LayerX1" in inputstring or "LayerU" in inputstring:
        return 0
    elif "LayerV" in inputstring or "LayerX2" in inputstring:
        return 1
    else:
        return -1

def module_to_int(inputstring):
    if "CSide" in inputstring:
        return 0
    elif "ASide" in inputstring:
        return 1
    else:
        return -1

def scale_lightness(rgb, scale_l):
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

colourname = 'technicolour'

label = "Tz"
limit = 0.1

grouping = ["CSide", "ASide"]
absolute=True

stations=["FT/T1","FT/T2","FT/T3"]
layers=["LayerX1","LayerU","LayerV","LayerX2"]

fig, axes = plt.subplots(nrows=len(stations),ncols=len(layers), figsize=(5*len(layers)+1, 3*len(stations)))
fig.tight_layout()

plt.figure()

def plot_RMSerrorsize_bygroup(align_output,plotted_alignables,title,filetitle,label,grouping,groupname,
                     absolute=True,fig=None,color="C0",nudge=0,iternum=-1):
    if fig is None:
        fig=plt.figure(figsize=(6,4))
    else:
        plt.figure(fig.number)
    xlist=[]
    for ii,group in enumerate(grouping):

        plt.title(title)
        if absolute:
            plt.ylabel(f'RMS error in abs. {label} (mm)' if "T" in label else f'abs. {label} (mrad)' if "R" in label else f'abs. {label}')
        else:
            plt.ylabel(f'RMS error in {label} (um)' if "T" in label else f'{label} (mrad)' if "R" in label else f'{label}')

        aggregate={"mean":0.0,"RMS":0.0}
        selected_alignables=[alignable for alignable in plotted_alignables if group in alignable]
        collected_list=[]
        for alignable in selected_alignables:
            if absolute:
                collected_list.append(abs(align_output[alignable][label][iternum]))
            else:
                collected_list.append(align_output[alignable][label][iternum])
        aggregate["mean"]=np.mean(collected_list)
        aggregate["RMS"]=np.sqrt(np.mean([(a-aggregate["mean"])**2 for a in collected_list]))
        #xlist.append(align_output[selected_alignables[0]]['x_global'][0])

        plt.scatter((align_output[selected_alignables[0]]['x_global'][-1]+nudge)/1000,1000*aggregate["RMS"],color=color,label=groupname if ii==0 else "_nolabel")
        plt.axhline(0,color='k',ls=':')

    return fig

def open_alignment(thisfile):
    with open(thisfile) as f:
        align_output=json.load(f)

    convergences=align_output.pop("converged")
    align_output.pop("total_chi2_vals", None)
    align_output.pop("total_chi2_dofs", None)
    #fix floats
    for alignable in align_output.keys():
        for label in labels+positions:
            align_output[alignable][label]=[float(ele.strip(',')) for ele in align_output[alignable][label]]
    return align_output

def choose_modules(align_output):
    plotted_alignables = align_output.keys()
    plotted_alignables = [alignable for alignable in plotted_alignables if "Mat" not in alignable]
#    plotted_alignables = [alignable for alignable in plotted_alignables if "Module" in alignable]
    plotted_alignables = [alignable for alignable in plotted_alignables if "Side" in alignable]
    return plotted_alignables

def plot_RMS_bygroup(align_output,plotted_alignables,title,filetitle,label,grouping,groupname,
                     absolute=True,fig=None,color="C0",nudge=0,iternum=-1,xaxis='z_global'):
    if fig is None:
        fig=plt.figure(figsize=(6,4))
    else:
        plt.figure(fig.number)
    xlist=[]
    for ii,group in enumerate(grouping):

        plt.title(title)
        if absolute:
            plt.ylabel(f'abs. {label} (mm)' if "T" in label else f'abs. {label} (mrad)' if "R" in label else f'abs. {label}')
        else:
            plt.ylabel(f'{label} (um)' if "T" in label else f'{label} (mrad)' if "R" in label else f'{label}')

        aggregate={"mean":0.0,"RMS":0.0}
        selected_alignables=[alignable for alignable in plotted_alignables if group in alignable]
        collected_list=[]
        for alignable in selected_alignables:
            if absolute:
                collected_list.append(abs(align_output[alignable][label][iternum]))
            else:
                collected_list.append(align_output[alignable][label][iternum])
        if len(collected_list) != 0:
            if isinstance(collected_list[0], str) == True:
                for st in range(len(collected_list)):
                    collected_list[st] = float(collected_list[st].strip(','))

            aggregate["mean"]=np.mean(collected_list)
            aggregate["RMS"]=np.sqrt(np.mean([(a-aggregate["mean"])**2 for a in collected_list]))
            plt.errorbar((align_output[selected_alignables[0]][xaxis][0]+nudge)/1000,1000*aggregate["mean"],yerr=1000*aggregate["RMS"],fmt=f"o{color}",label=groupname if ii==0 else "_nolabel")
            plt.axhline(0,color='k',ls=':')
        else:
            continue
    return fig

################### 2019 and 2020 data #####################
thisfile = "../../../output_alignment/20211110/MU/2019_data/AlignmentResults/parsedlog.json"
align_out_2019_5000MU = open_alignment(thisfile)
plotted_2019_5000MU = choose_modules(align_out_2019_5000MU)

thisfile = "../../../output_alignment/20211110/MU/2020_data/AlignmentResults/parsedlog.json"
align_out_2020_5000MU = open_alignment(thisfile)
plotted_2020_5000MU = choose_modules(align_out_2020_5000MU)

thisfile = "../../../output_alignment/20211110/MD/2019_data/AlignmentResults/parsedlog.json"
align_out_2019_5000MD = open_alignment(thisfile)
plotted_2019_5000MD = choose_modules(align_out_2019_5000MD)

thisfile = "../../../output_alignment/20211110/MD/2020_data/AlignmentResults/parsedlog.json"
align_out_2020_5000MD = open_alignment(thisfile)
plotted_2020_5000MD = choose_modules(align_out_2020_5000MD)

########################## needed data: ###########################
# n10 flo data
thisfile = "../output/data/Flo10/AlignmentResults/parsedlog.json"
align_out_Flo10_5000MU = open_alignment(thisfile)
plotted_Flo10_5000MU = choose_modules(align_out_Flo10_5000MU)

# config5, n10, e5000
thisfile = "../output/data/config5/backup/AlignmentResults/parsedlog.json"
align_out_config5_5000MU = open_alignment(thisfile)
plotted_config5_5000MU = choose_modules(align_out_config5_5000MU)

# config5, n10, e10000
thisfile = "../output/data/config5_10k/AlignmentResults/parsedlog.json"
align_out_config5_10000MU = open_alignment(thisfile)
plotted_config5_10000MU = choose_modules(align_out_config5_10000MU)

# config5, magdowwn, 2019
thisfile = "../output/data/c5_magdowwn_2019/AlignmentResults/parsedlog.json"
align_out_config5_5000MD_2019 = open_alignment(thisfile)
plotted_config5_5000MD_2019 = choose_modules(align_out_config5_5000MD_2019)

# config5, magdown, 2020
thisfile = "../output/data/c5_magdowwn_2020/AlignmentResults/parsedlog.json"
align_out_config5_5000MD_2020 = open_alignment(thisfile)
plotted_config5_5000MD_2020 = choose_modules(align_out_config5_5000MD_2020)

# misaligned data
thisfile = "../output/misalignment_runs/0/AlignmentResults/parsedlog.json"
align_out_T50_5000MU = open_alignment(thisfile)
plotted_T50_5000MU = choose_modules(align_out_T50_5000MU)

thisfile = "../output/misalignment_runs/1/AlignmentResults/parsedlog.json"
align_out_T200_5000MU = open_alignment(thisfile)
plotted_T200_5000MU = choose_modules(align_out_T200_5000MU)

thisfile = "../output/misalignment_runs/2/AlignmentResults/parsedlog.json"
align_out_R50_5000MU = open_alignment(thisfile)
plotted_R50_5000MU = choose_modules(align_out_R50_5000MU)

thisfile = "../output/misalignment_runs/3/AlignmentResults/parsedlog.json"
align_out_R100_5000MU = open_alignment(thisfile)
plotted_R100_5000MU = choose_modules(align_out_R100_5000MU)

############## plot all the alignment plots from "1" to "10" #########################
thisfile = "../output/alignment_runs/done_runs/HE_runs/0/AlignmentResults/parsedlog.json"
align_out_done0 = open_alignment(thisfile)
plotted_done0 = choose_modules(align_out_done0)

thisfile = "../output/alignment_runs/done_runs/HE_runs/1/AlignmentResults/parsedlog.json"
align_out_done1 = open_alignment(thisfile)
plotted_done1 = choose_modules(align_out_done1)

thisfile = "../output/alignment_runs/done_runs/HE_runs/2/AlignmentResults/parsedlog.json"
align_out_done2 = open_alignment(thisfile)
plotted_done2 = choose_modules(align_out_done2)

thisfile = "../output/alignment_runs/done_runs/HE_runs/3/AlignmentResults/parsedlog.json"
align_out_done3 = open_alignment(thisfile)
plotted_done3 = choose_modules(align_out_done3)

thisfile = "../output/alignment_runs/done_runs/HE_runs/4/AlignmentResults/parsedlog.json"
align_out_done4 = open_alignment(thisfile)
plotted_done4 = choose_modules(align_out_done4)

thisfile = "../output/alignment_runs/done_runs/HE_runs/5/AlignmentResults/parsedlog.json"
align_out_done5 = open_alignment(thisfile)
plotted_done5 = choose_modules(align_out_done5)

thisfile = "../output/alignment_runs/done_runs/HE_runs/6/AlignmentResults/parsedlog.json"
align_out_done6 = open_alignment(thisfile)
plotted_done6 = choose_modules(align_out_done6)

thisfile = "../output/alignment_runs/done_runs/HE_runs/7/AlignmentResults/parsedlog.json"
align_out_done7 = open_alignment(thisfile)
plotted_done7 = choose_modules(align_out_done7)

thisfile = "../output/alignment_runs/done_runs/HE_runs/8/AlignmentResults/parsedlog.json"
align_out_done8 = open_alignment(thisfile)
plotted_done8 = choose_modules(align_out_done8)

thisfile = "../output/alignment_runs/done_runs/HE_runs/9/AlignmentResults/parsedlog.json"
align_out_done9 = open_alignment(thisfile)
plotted_done9 = choose_modules(align_out_done9)

thisfile = "../output/alignment_runs/done_runs/HE_runs/10/AlignmentResults/parsedlog.json"
align_out_done10 = open_alignment(thisfile)
plotted_done10 = choose_modules(align_out_done10)

thisfile = "../output/alignment_runs/done_runs/HE_runs/11/AlignmentResults/parsedlog.json"
align_out_done11 = open_alignment(thisfile)
plotted_done11 = choose_modules(align_out_done11)

thisfile = "../output/alignment_runs/done_runs/HE_runs/12/iter20/AlignmentResults/parsedlog.json"
align_out_done12 = open_alignment(thisfile)
plotted_done12 = choose_modules(align_out_done12)

#####################
# config5_2019 again, 5000 events
thisfile = "../output/data/config5_2019/AlignmentResults/parsedlog.json"
align_out_config5_2019 = open_alignment(thisfile)
plotted_config5_2019 = choose_modules(align_out_config5_2019)

# config5_2020_fulliter, 10000 events


############## END DATA #################

grouping=["LayerX1","LayerU","LayerV","LayerX2"]
newgrouping=[]
for station in ["T1","T2","T3"]:
    for group in grouping:
        newgrouping.append(station+group)

filetitle = "RMSLayers"

###############################################################
outdir = "donefiles"

#12
plotted_alignables=[alignable for alignable in plotted_done12 if "CSide" in alignable]
fig81=plot_RMS_bygroup(align_out_done12, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"CSide of layers, withModules, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done12 if "ASide" in alignable]
fig81=plot_RMS_bygroup(align_out_done12, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=19,fig=fig81,color="C2",nudge=20)
plt.ylabel(r"Mean x translation ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Tx_{fig81.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done12 if "CSide" in alignable]
fig82=plot_RMS_bygroup(align_out_done12, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"CSide of layers, withModules, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done12 if "ASide" in alignable]
fig82=plot_RMS_bygroup(align_out_done12, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=19,fig=fig82,color="C2",nudge=20)
plt.ylabel(r"Mean z translation ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Tz_{fig82.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done12 if "CSide" in alignable]
fig83=plot_RMS_bygroup(align_out_done12, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"CSide of layers, withModules, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done12 if "ASide" in alignable]
fig83=plot_RMS_bygroup(align_out_done12, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=19,fig=fig83,color="C2",nudge=20)
plt.ylabel(r"Mean x rotation (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Rx_{fig83.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done12 if "CSide" in alignable]
fig84=plot_RMS_bygroup(align_out_done12, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"CSide of layers, withModules, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done12 if "ASide" in alignable]
fig84=plot_RMS_bygroup(align_out_done12, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=19,fig=fig84,color="C2",nudge=20)
plt.ylabel(r"Mean z rotation (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Rz_{fig84.number}.pdf")

# 0
fig17=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig17=plot_RMS_bygroup(align_out_done0, plotted_done0,title_t,filetitle,"Tx",newgrouping,"HE, config5, 2020, TxTzRxRz",absolute=False,iternum=9,fig=fig17,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done0_Tx_5000MU_{fig17.number}.pdf")

fig18=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig18=plot_RMS_bygroup(align_out_done0, plotted_done0,title_t,filetitle,"Tz",newgrouping,"HE, config5, 2020, TxTzRxRz",absolute=False,iternum=9,fig=fig18,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done0_Tz_5000MU_{fig18.number}.pdf")

fig19=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig19=plot_RMS_bygroup(align_out_done0, plotted_done0,title_r,filetitle,"Rx",newgrouping,"HE, config5, 2020, TxTzRxRz",absolute=False,iternum=9,fig=fig19,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done0_Rx_{fig19.number}.pdf")

fig20=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig20=plot_RMS_bygroup(align_out_done0, plotted_done0,title_r,filetitle,"Rz",newgrouping,"HE, config5, 2020, TxTzRxRz",absolute=False,iternum=9,fig=fig20,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done0_Rz_{fig20.number}.pdf")

# 1
fig21=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig21=plot_RMS_bygroup(align_out_done1, plotted_done1,title_t,filetitle,"Tx",newgrouping,"HE, config5, 2020, TxRz",absolute=False,iternum=9,fig=fig21,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done1_Tx_{fig21.number}.pdf")

fig22=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig22=plot_RMS_bygroup(align_out_done1, plotted_done1,title_t,filetitle,"Tz",newgrouping,"HE, config5, 2020, TxRz",absolute=False,iternum=9,fig=fig22,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done1_Tz_{fig22.number}.pdf")

fig23=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig23=plot_RMS_bygroup(align_out_done1, plotted_done1,title_r,filetitle,"Rx",newgrouping,"HE, config5, 2020, TxRz",absolute=False,iternum=9,fig=fig23,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done1_Rx_{fig23.number}.pdf")

fig24=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig24=plot_RMS_bygroup(align_out_done1, plotted_done1,title_r,filetitle,"Rz",newgrouping,"HE, config5, 2020, TxRz",absolute=False,iternum=9,fig=fig24,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done1_Rz_{fig24.number}.pdf")

# 2
fig25=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig25=plot_RMS_bygroup(align_out_done2, plotted_done2,title_t,filetitle,"Tx",newgrouping,"HE, c5, ModulesAligned, 2020, TxTzRxRz",absolute=False,iternum=9,fig=fig25,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done2_Tx_{fig25.number}.pdf")

fig26=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig26=plot_RMS_bygroup(align_out_done2, plotted_done2,title_t,filetitle,"Tz",newgrouping,"HE, c5, ModulesAligned, 2020, TxTzRxRz",absolute=False,iternum=9,fig=fig26,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done2_Tz_{fig26.number}.pdf")

fig27=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig27=plot_RMS_bygroup(align_out_done2, plotted_done2,title_r,filetitle,"Rx",newgrouping,"HE, c5, ModulesAligned, 2020, TxTzRxRz",absolute=False,iternum=9,fig=fig27,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done2_Rx_{fig27.number}.pdf")

fig28=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig28=plot_RMS_bygroup(align_out_done2, plotted_done2,title_r,filetitle,"Rz",newgrouping,"HE, c5, ModulesAligned, 2020, TxTzRxRz",absolute=False,iternum=9,fig=fig28,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done2_Rz_{fig28.number}.pdf")

# 3
fig29=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig29=plot_RMS_bygroup(align_out_done3, plotted_done3,title_t,filetitle,"Tx",newgrouping,"HE, c5, ModulesAligned, 2020, TxRz",absolute=False,iternum=9,fig=fig29,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done3_Tx_{fig29.number}.pdf")

fig30=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig30=plot_RMS_bygroup(align_out_done3, plotted_done3,title_t,filetitle,"Tz",newgrouping,"HE, c5, ModulesAligned, 2020, TxRz",absolute=False,iternum=9,fig=fig30,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done3_Tz_{fig30.number}.pdf")

fig31=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig31=plot_RMS_bygroup(align_out_done3, plotted_done3,title_r,filetitle,"Rx",newgrouping,"HE, c5, ModulesAligned, 2020, TxRz",absolute=False,iternum=9,fig=fig31,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done3_Rx_{fig31.number}.pdf")

fig32=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig32=plot_RMS_bygroup(align_out_done3, plotted_done3,title_r,filetitle,"Rz",newgrouping,"HE, c5, ModulesAligned, 2020, TxRz",absolute=False,iternum=9,fig=fig32,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done3_Rz_{fig32.number}.pdf")

# 4
fig33=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig33=plot_RMS_bygroup(align_out_done4, plotted_done4,title_t,filetitle,"Tx",newgrouping,"HE, config5, 2020, TxTzRxRz, 16 iter",absolute=False,iternum=9,fig=fig33,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done4_Tx_{fig33.number}.pdf")

fig34=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig34=plot_RMS_bygroup(align_out_done4, plotted_done4,title_t,filetitle,"Tz",newgrouping,"HE, config5, 2020, TxTzRxRz, 16 iter",absolute=False,iternum=9,fig=fig34,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done4_Tz_{fig34.number}.pdf")

fig35=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig35=plot_RMS_bygroup(align_out_done4, plotted_done4,title_r,filetitle,"Rx",newgrouping,"HE, config5, 2020, TxTzRxRz, 16 iter",absolute=False,iternum=9,fig=fig35,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done4_Rx_{fig35.number}.pdf")

fig36=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig36=plot_RMS_bygroup(align_out_done4, plotted_done4,title_r,filetitle,"Rz",newgrouping,"HE, config5, 2020, TxTzRxRz, 16 iter",absolute=False,iternum=9,fig=fig36,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done4_Rz_{fig36.number}.pdf")

# 5
fig45=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig45=plot_RMS_bygroup(align_out_done7, plotted_done7,title_t,filetitle,"Tx",newgrouping,"HE, config5, 2019, 6500, TxRz",absolute=False,iternum=9,fig=fig45,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done7_Tx_{fig45.number}.pdf")

fig46=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig46=plot_RMS_bygroup(align_out_done7, plotted_done7,title_t,filetitle,"Tz",newgrouping,"HE, config5, 2019, 6500, TxRz",absolute=False,iternum=9,fig=fig46,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done7_Tz_{fig46.number}.pdf")

fig47=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig47=plot_RMS_bygroup(align_out_done7, plotted_done7,title_r,filetitle,"Rx",newgrouping,"HE, config5, 2019, 6500, TxRz",absolute=False,iternum=9,fig=fig47,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done7_Rx_{fig47.number}.pdf")

fig48=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig48=plot_RMS_bygroup(align_out_done7, plotted_done7,title_r,filetitle,"Rz",newgrouping,"HE, config5, 2019, 6500, TxRz",absolute=False,iternum=9,fig=fig48,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done7_Rz_{fig48.number}.pdf")

# 6
fig49=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig49=plot_RMS_bygroup(align_out_done8, plotted_done8,title_t,filetitle,"Tx",newgrouping,"HE, config5 with Modules, 2019, 6500, TxRz",absolute=False,iternum=9,fig=fig49,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done8_Tx_{fig49.number}.pdf")

fig50=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig50=plot_RMS_bygroup(align_out_done8, plotted_done8,title_t,filetitle,"Tz",newgrouping,"HE, config5 with Modules, 2019, 6500, TxRz",absolute=False,iternum=9,fig=fig50,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done8_Tz_{fig50.number}.pdf")

fig51=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig51=plot_RMS_bygroup(align_out_done8, plotted_done8,title_r,filetitle,"Rx",newgrouping,"HE, config5 with Modules, 2019, 6500, TxRz",absolute=False,iternum=9,fig=fig51,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done8_Rx_{fig51.number}.pdf")

fig52=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig52=plot_RMS_bygroup(align_out_done8, plotted_done8,title_r,filetitle,"Rz",newgrouping,"HE, config5 with Modules, 2019, 6500, TxRz",absolute=False,iternum=9,fig=fig52,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done8_Rz_{fig52.number}.pdf")

# 7
fig37=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig37=plot_RMS_bygroup(align_out_done5, plotted_done5,title_t,filetitle,"Tx",newgrouping,"HE, config5, 2019, 30k, TxRz",absolute=False,iternum=9,fig=fig37,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done5_Tx_{fig37.number}.pdf")

fig38=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig38=plot_RMS_bygroup(align_out_done5, plotted_done5,title_t,filetitle,"Tz",newgrouping,"HE, config5, 2019, 30k, TxRz",absolute=False,iternum=9,fig=fig38,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done5_Tz_{fig38.number}.pdf")

fig39=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig39=plot_RMS_bygroup(align_out_done5, plotted_done5,title_r,filetitle,"Rx",newgrouping,"HE, config5, 2019, 30k, TxRz",absolute=False,iternum=9,fig=fig39,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done5_Rx_{fig39.number}.pdf")

fig40=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig40=plot_RMS_bygroup(align_out_done5, plotted_done5,title_r,filetitle,"Rz",newgrouping,"HE, config5, 2019, 30k, TxRz",absolute=False,iternum=9,fig=fig40,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done5_Rz_{fig40.number}.pdf")

#8
fig41=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig41=plot_RMS_bygroup(align_out_done6, plotted_done6,title_t,filetitle,"Tx",newgrouping,"HE, config5, 2019, 30k, TxRz, ModulesAligned",absolute=False,iternum=9,fig=fig41,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done6_Tx_{fig41.number}.pdf")

fig42=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig42=plot_RMS_bygroup(align_out_done6, plotted_done6,title_t,filetitle,"Tz",newgrouping,"HE, config5, 2019, 30k, TxRz, ModulesAligned",absolute=False,iternum=9,fig=fig42,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done6_Tz_{fig42.number}.pdf")

fig43=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig43=plot_RMS_bygroup(align_out_done6, plotted_done6,title_r,filetitle,"Rx",newgrouping,"HE, config5, 2019, 30k, TxRz, ModulesAligned",absolute=False,iternum=9,fig=fig43,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done6_Rx_{fig43.number}.pdf")

fig44=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig44=plot_RMS_bygroup(align_out_done6, plotted_done6,title_r,filetitle,"Rz",newgrouping,"HE, config5, 2019, 30k, TxRz, ModulesAligned",absolute=False,iternum=9,fig=fig44,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done6_Rz_{fig44.number}.pdf")

#9
fig53=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig53=plot_RMS_bygroup(align_out_done9, plotted_done9,title_t,filetitle,"Tx",newgrouping,"HE, config5, 2020, TxRz, 10k events",absolute=False,iternum=9,fig=fig53,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done0_Tx_5000MU_{fig53.number}.pdf")

fig54=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig54=plot_RMS_bygroup(align_out_done9, plotted_done9,title_t,filetitle,"Tz",newgrouping,"HE, config5, 2020, TxRz, 10k events",absolute=False,iternum=9,fig=fig54,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done0_Tz_5000MU_{fig54.number}.pdf")

fig55=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig55=plot_RMS_bygroup(align_out_done9, plotted_done9,title_r,filetitle,"Rx",newgrouping,"HE, config5, 2020, TxRz, 10k events",absolute=False,iternum=9,fig=fig55,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done0_Rx_{fig55.number}.pdf")

fig56=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig56=plot_RMS_bygroup(align_out_done9, plotted_done9,title_r,filetitle,"Rz",newgrouping,"HE, config5, 2020, TxRz, 10k events",absolute=False,iternum=9,fig=fig56,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done0_Rz_{fig56.number}.pdf")

#10
fig57=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig57=plot_RMS_bygroup(align_out_done10, plotted_done10,title_t,filetitle,"Tx",newgrouping,"HE, config5, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9,fig=fig57,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done1_Tx_{fig57.number}.pdf")

fig58=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig58=plot_RMS_bygroup(align_out_done10, plotted_done10,title_t,filetitle,"Tz",newgrouping,"HE, config5, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9,fig=fig58,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done1_Tz_{fig58.number}.pdf")

fig59=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig59=plot_RMS_bygroup(align_out_done10, plotted_done10,title_r,filetitle,"Rx",newgrouping,"HE, config5, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9,fig=fig59,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done1_Rx_{fig59.number}.pdf")

fig60=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig60=plot_RMS_bygroup(align_out_done10, plotted_done10,title_r,filetitle,"Rz",newgrouping,"HE, config5, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9,fig=fig60,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done1_Rz_{fig60.number}.pdf")

#11
fig61=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig61=plot_RMS_bygroup(align_out_done11, plotted_done11,title_t,filetitle,"Tx",newgrouping,"HE, c5, withModules, 2020, TxRz, 10k, events",absolute=False,iternum=9,fig=fig61,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done2_Tx_{fig61.number}.pdf")

fig62=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig62=plot_RMS_bygroup(align_out_done11, plotted_done11,title_t,filetitle,"Tz",newgrouping,"HE, c5, withModules, 2020, TxRz, 10k events",absolute=False,iternum=9,fig=fig62,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done2_Tz_{fig62.number}.pdf")

fig63=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig63=plot_RMS_bygroup(align_out_done11, plotted_done11,title_r,filetitle,"Rx",newgrouping,"HE, c5, withModules, 2020, TxRz, 10k events",absolute=False,iternum=9,fig=fig63,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done2_Rx_{fig63.number}.pdf")

fig64=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig64=plot_RMS_bygroup(align_out_done11, plotted_done11,title_r,filetitle,"Rz",newgrouping,"HE, c5, withModules, 2020, TxRz, 10k events",absolute=False,iternum=9,fig=fig64,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done2_Rz_{fig64.number}.pdf")

#12
fig65=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig65=plot_RMS_bygroup(align_out_done12, plotted_done12,title_t,filetitle,"Tx",newgrouping,"HE, c5, withModules, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9,fig=fig65,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done3_Tx_{fig65.number}.pdf")

fig66=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig66=plot_RMS_bygroup(align_out_done12, plotted_done12,title_t,filetitle,"Tz",newgrouping,"HE, c5, withModules, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9,fig=fig66,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done3_Tz_{fig66.number}.pdf")

fig67=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig67=plot_RMS_bygroup(align_out_done12, plotted_done12,title_r,filetitle,"Rx",newgrouping,"HE, c5, withModules, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9,fig=fig67,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done3_Rx_{fig67.number}.pdf")

fig68=plot_RMS_bygroup(align_out_config5_10000MU, plotted_config5_10000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 10000MU",absolute=False,iternum=9)
fig68=plot_RMS_bygroup(align_out_done12, plotted_done12,title_r,filetitle,"Rz",newgrouping,"HE, c5, withModules, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9,fig=fig68,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done3_Rz_{fig68.number}.pdf")

# C-A side plots for 9 to 12
#9
plotted_alignables=[alignable for alignable in plotted_done9 if "CSide" in alignable]
fig69=plot_RMS_bygroup(align_out_done9, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"CSide of layers, 2020, TxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done9 if "ASide" in alignable]
fig69=plot_RMS_bygroup(align_out_done9, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig69,color="C2",nudge=20)
plt.ylabel(r"Mean x translation ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Tx_{fig69.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done9 if "CSide" in alignable]
fig70=plot_RMS_bygroup(align_out_done9, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"CSide of layers, 2020, TxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done9 if "ASide" in alignable]
fig70=plot_RMS_bygroup(align_out_done9, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig70,color="C2",nudge=20)
plt.ylabel(r"Mean z translation ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Tz_{fig70.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done9 if "CSide" in alignable]
fig71=plot_RMS_bygroup(align_out_done9, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"CSide of layers, 2020, TxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done9 if "ASide" in alignable]
fig71=plot_RMS_bygroup(align_out_done9, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig71,color="C2",nudge=20)
plt.ylabel(r"Mean x rotation (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Rx_{fig71.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done9 if "CSide" in alignable]
fig72=plot_RMS_bygroup(align_out_done9, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"CSide of layers, 2020, TxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done9 if "ASide" in alignable]
fig72=plot_RMS_bygroup(align_out_done9, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig72,color="C2",nudge=20)
plt.ylabel(r"Mean z rotation (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Rz_{fig72.number}.pdf")

#10
plotted_alignables=[alignable for alignable in plotted_done10 if "CSide" in alignable]
fig73=plot_RMS_bygroup(align_out_done10, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"CSide of layers, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done10 if "ASide" in alignable]
fig73=plot_RMS_bygroup(align_out_done10, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig73,color="C2",nudge=20)
plt.ylabel(r"Mean x translation ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Tx_{fig73.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done10 if "CSide" in alignable]
fig74=plot_RMS_bygroup(align_out_done10, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"CSide of layers, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done10 if "ASide" in alignable]
fig74=plot_RMS_bygroup(align_out_done10, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig74,color="C2",nudge=20)
plt.ylabel(r"Mean z translation ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Tz_{fig74.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done10 if "CSide" in alignable]
fig75=plot_RMS_bygroup(align_out_done10, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"CSide of layers, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done10 if "ASide" in alignable]
fig75=plot_RMS_bygroup(align_out_done10, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig75,color="C2",nudge=20)
plt.ylabel(r"Mean x rotation (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Rx_{fig75.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done10 if "CSide" in alignable]
fig76=plot_RMS_bygroup(align_out_done10, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"CSide of layers, 2020, TxTzRxRz, 10k events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done10 if "ASide" in alignable]
fig76=plot_RMS_bygroup(align_out_done10, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig76,color="C2",nudge=20)
plt.ylabel(r"Mean z rotation (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Rz_{fig76.number}.pdf")

#11
plotted_alignables=[alignable for alignable in plotted_done11 if "CSide" in alignable]
fig77=plot_RMS_bygroup(align_out_done11, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"CSide of layers, withModules, 2020, TxRz, 10k, events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done11 if "ASide" in alignable]
fig77=plot_RMS_bygroup(align_out_done11, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig77,color="C2",nudge=20)
plt.ylabel(r"Mean x translation ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Tx_{fig77.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done11 if "CSide" in alignable]
fig78=plot_RMS_bygroup(align_out_done11, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"CSide of layers, withModules, 2020, TxRz, 10k, events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done11 if "ASide" in alignable]
fig78=plot_RMS_bygroup(align_out_done11, plotted_alignables,title_t,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig78,color="C2",nudge=20)
plt.ylabel(r"Mean z translation ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Tz_{fig78.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done11 if "CSide" in alignable]
fig79=plot_RMS_bygroup(align_out_done11, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"CSide of layers, withModules, 2020, TxRz, 10k, events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done11 if "ASide" in alignable]
fig79=plot_RMS_bygroup(align_out_done11, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig79,color="C2",nudge=20)
plt.ylabel(r"Mean x rotation (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Rx_{fig79.number}.pdf")

plotted_alignables=[alignable for alignable in plotted_done11 if "CSide" in alignable]
fig80=plot_RMS_bygroup(align_out_done11, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"CSide of layers, withModules, 2020, TxRz, 10k, events",absolute=False,iternum=9)
plotted_alignables=[alignable for alignable in plotted_done11 if "ASide" in alignable]
fig80=plot_RMS_bygroup(align_out_done11, plotted_alignables,title_r,filetitle,"Tx",newgrouping,"ASide of layers",absolute=False,iternum=9,fig=fig80,color="C2",nudge=20)
plt.ylabel(r"Mean z rotation (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/done9_CA_Rz_{fig80.number}.pdf")

######### test earlier config5 with now config5_2019
# 85-88
fig85=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2019 data, 5000MU",absolute=False,iternum=9)
fig85=plot_RMS_bygroup(align_out_config5_2019, plotted_config5_2019,title_t,filetitle,"Tx",newgrouping,"config5, 2019, 5000MU, retest",absolute=False,iternum=9,fig=fig85,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/c5_retest_Tx_{fig85.number}.pdf")

fig86=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2019 data, 5000MU",absolute=False,iternum=9)
fig86=plot_RMS_bygroup(align_out_config5_2019, plotted_config5_2019,title_t,filetitle,"Tz",newgrouping,"config5, 2019, retest",absolute=False,iternum=9,fig=fig86,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/c5_retest_Tz_{fig86.number}.pdf")

fig87=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2019 data, 5000MU",absolute=False,iternum=9)
fig87=plot_RMS_bygroup(align_out_config5_2019, plotted_config5_2019,title_r,filetitle,"Rx",newgrouping,"config5, 2019, retest",absolute=False,iternum=9,fig=fig87,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/c5_retest_Rx_{fig87.number}.pdf")

fig88=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2019 data, 5000MU",absolute=False,iternum=9)
fig88=plot_RMS_bygroup(align_out_config5_2019, plotted_config5_2019,title_r,filetitle,"Rz",newgrouping,"config5, 2019, retest",absolute=False,iternum=9,fig=fig88,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/c5_retest_Rz_{fig88.number}.pdf")


plt.show()
