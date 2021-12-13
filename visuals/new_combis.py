import matplotlib.pyplot as plt

#for recolor
import matplotlib.colors as mplcolor
import colorsys as colorsys

import json
import numpy as np

outdir = "outfiles"

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

    #fix floats
    for alignable in align_output.keys():
        for label in labels+positions:
            align_output[alignable][label]=[float(ele.strip(',')) for ele in align_output[alignable][label]]
    return align_output

def choose_modules(align_output):
    plotted_alignables = align_output.keys()
    plotted_alignables = [alignable for alignable in plotted_alignables if "Mat" not in alignable]
    #plotted_alignables = [alignable for alignable in plotted_alignables if "Module" in alignable]
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
# converts string numbers to real floats
#        print(collected_list)
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

################ DATA #################
# florian5000, n = 12
#thisfile = "../output/data/AlignmentResults/parsedlog.json"
#align_out_Flo12_5000MU = open_alignment(thisfile)
#plotted_Flo12_5000MU = choose_modules(align_out_Flo12_5000MU)

#
#thisfile = "../output/alignment_runs/0/AlignmentResults/parsedlog.json"
#align_out_conf5_Rz_5000MU = open_alignment(thisfile)
#plotted_conf5_Rz_5000MU = choose_modules(align_out_conf5_Rz_5000MU)

#
#thisfile = "../output/alignment_runs/1/AlignmentResults/parsedlog.json"
#align_out_noTx_5000MU = open_alignment(thisfile)
#plotted_noTx_5000MU = choose_modules(align_out_noTx_5000MU)

#
#thisfile = "../output/alignment_runs/2/AlignmentResults/parsedlog.json"
#align_out_TxRz_5000MU = open_alignment(thisfile)
#plotted_TxRz_5000MU = choose_modules(align_out_TxRz_5000MU)

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
thisfile = "../output/data/config5/AlignmentResults/parsedlog.json"
align_out_config5_5000MU = open_alignment(thisfile)
plotted_config5_5000MU = choose_modules(align_out_config5_5000MU)

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

# 8 misaligned translations

thisfile = "../output/misalignment_runs/0/AlignmentResults/parsedlog.json"
align_out_T0_5000MU = open_alignment(thisfile)
plotted_T0_5000MU = choose_modules(align_out_T0_5000MU)

thisfile = "../output/misalignment_runs/1/AlignmentResults/parsedlog.json"
align_out_T1_5000MU = open_alignment(thisfile)
plotted_T1_5000MU = choose_modules(align_out_T1_5000MU)

thisfile = "../output/misalignment_runs/2/AlignmentResults/parsedlog.json"
align_out_T2_5000MU = open_alignment(thisfile)
plotted_T2_5000MU = choose_modules(align_out_T2_5000MU)

thisfile = "../output/misalignment_runs/3/AlignmentResults/parsedlog.json"
align_out_T3_5000MU = open_alignment(thisfile)
plotted_T3_5000MU = choose_modules(align_out_T3_5000MU)

thisfile = "../output/misalignment_runs/4/AlignmentResults/parsedlog.json"
align_out_T4_5000MU = open_alignment(thisfile)
plotted_T4_5000MU = choose_modules(align_out_T4_5000MU)

thisfile = "../output/misalignment_runs/5/AlignmentResults/parsedlog.json"
align_out_T5_5000MU = open_alignment(thisfile)
plotted_T5_5000MU = choose_modules(align_out_T5_5000MU)

thisfile = "../output/misalignment_runs/6/AlignmentResults/parsedlog.json"
align_out_T6_5000MU = open_alignment(thisfile)
plotted_T6_5000MU = choose_modules(align_out_T6_5000MU)

thisfile = "../output/misalignment_runs/7/AlignmentResults/parsedlog.json"
align_out_T7_5000MU = open_alignment(thisfile)
plotted_T7_5000MU = choose_modules(align_out_T7_5000MU)

### misaligned high momentum tracks
thisfile = "../output/misalignment_runs/done_runs/HETTracks_8x/0/AlignmentResults/parsedlog.json"
align_out_HT0_5000MU = open_alignment(thisfile)
plotted_HT0_5000MU = choose_modules(align_out_HT0_5000MU)

thisfile = "../output/misalignment_runs/done_runs/HETTracks_8x/1/AlignmentResults/parsedlog.json"
align_out_HT1_5000MU = open_alignment(thisfile)
plotted_HT1_5000MU = choose_modules(align_out_HT1_5000MU)

thisfile = "../output/misalignment_runs/done_runs/HETTracks_8x/2/AlignmentResults/parsedlog.json"
align_out_HT2_5000MU = open_alignment(thisfile)
plotted_HT2_5000MU = choose_modules(align_out_HT2_5000MU)

thisfile = "../output/misalignment_runs/done_runs/HETTracks_8x/3/AlignmentResults/parsedlog.json"
align_out_HT3_5000MU = open_alignment(thisfile)
plotted_HT3_5000MU = choose_modules(align_out_HT3_5000MU)

thisfile = "../output/misalignment_runs/done_runs/HETTracks_8x/4/AlignmentResults/parsedlog.json"
align_out_HT4_5000MU = open_alignment(thisfile)
plotted_HT4_5000MU = choose_modules(align_out_HT4_5000MU)

thisfile = "../output/misalignment_runs/done_runs/HETTracks_8x/5/AlignmentResults/parsedlog.json"
align_out_HT5_5000MU = open_alignment(thisfile)
plotted_HT5_5000MU = choose_modules(align_out_HT5_5000MU)

thisfile = "../output/misalignment_runs/done_runs/HETTracks_8x/6/AlignmentResults/parsedlog.json"
align_out_HT6_5000MU = open_alignment(thisfile)
plotted_HT6_5000MU = choose_modules(align_out_HT6_5000MU)

thisfile = "../output/misalignment_runs/done_runs/HETTracks_8x/7/AlignmentResults/parsedlog.json"
align_out_HT7_5000MU = open_alignment(thisfile)
plotted_HT7_5000MU = choose_modules(align_out_HT7_5000MU)


############## END DATA #################

grouping=["LayerX1","LayerU","LayerV","LayerX2"]
newgrouping=[]
for station in ["T1","T2","T3"]:
    for group in grouping:
        newgrouping.append(station+group)

filetitle = "RMSLayers"

############## PLOTS #############
## misalignment runs
## T50
fig100=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_t,filetitle,"Tx",newgrouping,"fullL, Tx, flo, 5000MU",absolute=False,iternum=9)
fig100=plot_RMS_bygroup(align_out_T50_5000MU, plotted_T50_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 50 mu translation, Tx, 5000MU",absolute=False,iternum=9,fig=fig100,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig100.number}.pdf")

fig101=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_t,filetitle,"Tz",newgrouping,"fullL, Tz, flo, 5000MU",absolute=False,iternum=9)
fig101=plot_RMS_bygroup(align_out_T50_5000MU, plotted_T50_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 50 mu translation, Tz, 5000MU",absolute=False,iternum=9,fig=fig101,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig101.number}.pdf")

fig102=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_r,filetitle,"Rx",newgrouping,"fullL, Rx, flo, 5000MU",absolute=False,iternum=9)
fig102=plot_RMS_bygroup(align_out_T50_5000MU, plotted_T50_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 50 mu translation, Rx, 5000MU",absolute=False,iternum=9,fig=fig102,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig102.number}.pdf")

fig103=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_r,filetitle,"Rz",newgrouping,"fullL, Rz, flo, 5000MU",absolute=False,iternum=9)
fig103=plot_RMS_bygroup(align_out_T50_5000MU, plotted_T50_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 50 mu translation, Rz, 5000MU",absolute=False,iternum=9,fig=fig103,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig103.number}.pdf")

## T200
fig104=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_t,filetitle,"Tx",newgrouping,"fullL, Tx, flo, 5000MU",absolute=False,iternum=9)
fig104=plot_RMS_bygroup(align_out_T200_5000MU, plotted_T200_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 200 mu translation, Tx, 5000MU",absolute=False,iternum=9,fig=fig104,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig104.number}.pdf")

fig105=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_t,filetitle,"Tz",newgrouping,"fullL, Tz, flo, 5000MU",absolute=False,iternum=9)
fig105=plot_RMS_bygroup(align_out_T200_5000MU, plotted_T200_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 200 mu translation, Tz, 5000MU",absolute=False,iternum=9,fig=fig105,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig105.number}.pdf")

fig106=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_r,filetitle,"Rx",newgrouping,"fullL, Rx, flo, 5000MU",absolute=False,iternum=9)
fig106=plot_RMS_bygroup(align_out_T200_5000MU, plotted_T200_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 200 mu translation, Rx, 5000MU",absolute=False,iternum=9,fig=fig106,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig106.number}.pdf")

fig107=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_r,filetitle,"Rz",newgrouping,"fullL, Rz, flo, 5000MU",absolute=False,iternum=9)
fig107=plot_RMS_bygroup(align_out_T200_5000MU, plotted_T200_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 200 mu translation, Rz, 5000MU",absolute=False,iternum=9,fig=fig107,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig107.number}.pdf")

## R50
fig108=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_t,filetitle,"Tx",newgrouping,"fullL, Tx, flo, 5000MU",absolute=False,iternum=9)
fig108=plot_RMS_bygroup(align_out_R50_5000MU, plotted_R50_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 50 mu rotation, Tx, 5000MU",absolute=False,iternum=9,fig=fig108,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig108.number}.pdf")

fig109=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_t,filetitle,"Tz",newgrouping,"fullL, Tz, flo, 5000MU",absolute=False,iternum=9)
fig109=plot_RMS_bygroup(align_out_R50_5000MU, plotted_R50_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 50 mu rotation, Tz, 5000MU",absolute=False,iternum=9,fig=fig109,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig101.number}.pdf")

fig110=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_r,filetitle,"Rx",newgrouping,"fullL, Rx, flo, 5000MU",absolute=False,iternum=9)
fig110=plot_RMS_bygroup(align_out_R50_5000MU, plotted_R50_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 50 mu rotation, Rx, 5000MU",absolute=False,iternum=9,fig=fig110,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig110.number}.pdf")

fig111=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_r,filetitle,"Rz",newgrouping,"fullL, Rz, flo, 5000MU",absolute=False,iternum=9)
fig111=plot_RMS_bygroup(align_out_R50_5000MU, plotted_R50_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 50 mu rotation, Rz, 5000MU",absolute=False,iternum=9,fig=fig111,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig111.number}.pdf")

## R100
fig112=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_t,filetitle,"Tx",newgrouping,"fullL, Tx, flo, 5000MU",absolute=False,iternum=9)
fig112=plot_RMS_bygroup(align_out_R100_5000MU, plotted_R100_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 100 mu rotation, Tx, 5000MU",absolute=False,iternum=9,fig=fig112,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig112.number}.pdf")

fig113=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_t,filetitle,"Tz",newgrouping,"fullL, Tz, flo, 5000MU",absolute=False,iternum=9)
fig113=plot_RMS_bygroup(align_out_R100_5000MU, plotted_R100_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 100 mu rotation, Tz, 5000MU",absolute=False,iternum=9,fig=fig113,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig113.number}.pdf")

fig114=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_r,filetitle,"Rx",newgrouping,"fullL, Rx, flo, 5000MU",absolute=False,iternum=9)
fig114=plot_RMS_bygroup(align_out_R100_5000MU, plotted_R100_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 100 mu rotation, Rx, 5000MU",absolute=False,iternum=9,fig=fig114,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig114.number}.pdf")

fig115=plot_RMS_bygroup(align_out_Flo10_5000MU, plotted_Flo10_5000MU, title_r,filetitle,"Rz",newgrouping,"fullL, Rz, flo, 5000MU",absolute=False,iternum=9)
fig115=plot_RMS_bygroup(align_out_R100_5000MU, plotted_R100_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 100 mu rotation, Rz, 5000MU",absolute=False,iternum=9,fig=fig115,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig115.number}.pdf")

############### 8 misaligned translations for cross cheecking ############
# 0
fig120=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig120=plot_RMS_bygroup(align_out_T0_5000MU, plotted_T0_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 200 mu translation, Tx, 5000MU",absolute=False,iternum=9,fig=fig120,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig120.number}.pdf")

fig121=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig121=plot_RMS_bygroup(align_out_T0_5000MU, plotted_T0_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 200 mu translation, Tz, 5000MU",absolute=False,iternum=9,fig=fig121,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig121.number}.pdf")

fig122=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig122=plot_RMS_bygroup(align_out_T0_5000MU, plotted_T0_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 200 mu translation, Rx, 5000MU",absolute=False,iternum=9,fig=fig122,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig122.number}.pdf")

fig123=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig123=plot_RMS_bygroup(align_out_T0_5000MU, plotted_T0_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 200 mu translation, Rz, 5000MU",absolute=False,iternum=9,fig=fig123,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig123.number}.pdf")

# 1
fig120=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig120=plot_RMS_bygroup(align_out_T1_5000MU, plotted_T1_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 100 mu translation, Tx, 5000MU",absolute=False,iternum=9,fig=fig120,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig120.number}.pdf")

fig121=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig121=plot_RMS_bygroup(align_out_T1_5000MU, plotted_T1_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 100 mu translation, Tz, 5000MU",absolute=False,iternum=9,fig=fig121,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig121.number}.pdf")

fig122=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig122=plot_RMS_bygroup(align_out_T1_5000MU, plotted_T1_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 100 mu translation, Rx, 5000MU",absolute=False,iternum=9,fig=fig122,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig122.number}.pdf")

fig123=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig123=plot_RMS_bygroup(align_out_T1_5000MU, plotted_T1_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 100 mu translation, Rz, 5000MU",absolute=False,iternum=9,fig=fig123,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig123.number}.pdf")

# 2
fig124=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig124=plot_RMS_bygroup(align_out_T2_5000MU, plotted_T2_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 50 mu translation, Tx, 5000MU",absolute=False,iternum=9,fig=fig124,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig124.number}.pdf")

fig125=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig125=plot_RMS_bygroup(align_out_T2_5000MU, plotted_T2_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 50 mu translation, Tz, 5000MU",absolute=False,iternum=9,fig=fig125,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig125.number}.pdf")

fig126=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig126=plot_RMS_bygroup(align_out_T2_5000MU, plotted_T2_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 50 mu translation, Rx, 5000MU",absolute=False,iternum=9,fig=fig126,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig126.number}.pdf")

fig127=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig127=plot_RMS_bygroup(align_out_T2_5000MU, plotted_T2_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 50 mu translation, Rz, 5000MU",absolute=False,iternum=9,fig=fig127,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig127.number}.pdf")

# 3
fig128=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig128=plot_RMS_bygroup(align_out_T3_5000MU, plotted_T3_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 10 mu translation, Tx, 5000MU",absolute=False,iternum=9,fig=fig128,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig128.number}.pdf")

fig129=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig129=plot_RMS_bygroup(align_out_T3_5000MU, plotted_T3_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 10 mu translation, Tz, 5000MU",absolute=False,iternum=9,fig=fig129,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig129.number}.pdf")

fig130=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig130=plot_RMS_bygroup(align_out_T3_5000MU, plotted_T3_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 10 mu translation, Rx, 5000MU",absolute=False,iternum=9,fig=fig130,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig130.number}.pdf")

fig131=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig131=plot_RMS_bygroup(align_out_T3_5000MU, plotted_T3_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 10 mu translation, Rz, 5000MU",absolute=False,iternum=9,fig=fig131,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig131.number}.pdf")

# 4
fig132=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig132=plot_RMS_bygroup(align_out_T4_5000MU, plotted_T4_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 1 mu translation, Tx, 5000MU",absolute=False,iternum=9,fig=fig132,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig132.number}.pdf")

fig133=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig133=plot_RMS_bygroup(align_out_T4_5000MU, plotted_T4_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 1 mu translation, Tz, 5000MU",absolute=False,iternum=9,fig=fig133,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig133.number}.pdf")

fig134=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig134=plot_RMS_bygroup(align_out_T4_5000MU, plotted_T4_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 1 mu translation, Rx, 5000MU",absolute=False,iternum=9,fig=fig134,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig134.number}.pdf")

fig135=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig135=plot_RMS_bygroup(align_out_T4_5000MU, plotted_T4_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 1 mu translation, Rz, 5000MU",absolute=False,iternum=9,fig=fig135,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig135.number}.pdf")

# 5
fig136=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig136=plot_RMS_bygroup(align_out_T5_5000MU, plotted_T5_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 0.1 mu translation, Tx, 5000MU",absolute=False,iternum=9,fig=fig136,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig136.number}.pdf")

fig137=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig137=plot_RMS_bygroup(align_out_T5_5000MU, plotted_T5_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 0.1 mu translation, Tz, 5000MU",absolute=False,iternum=9,fig=fig137,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig137.number}.pdf")

fig138=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig138=plot_RMS_bygroup(align_out_T5_5000MU, plotted_T5_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 0.1 mu translation, Rx, 5000MU",absolute=False,iternum=9,fig=fig138,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig138.number}.pdf")

fig139=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig139=plot_RMS_bygroup(align_out_T5_5000MU, plotted_T5_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 0.1 mu translation, Rz, 5000MU",absolute=False,iternum=9,fig=fig139,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig139.number}.pdf")

# 6
fig140=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig140=plot_RMS_bygroup(align_out_T6_5000MU, plotted_T6_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 0.01 mu translation, Tx, 5000MU",absolute=False,iternum=9,fig=fig140,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig140.number}.pdf")

fig141=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig141=plot_RMS_bygroup(align_out_T6_5000MU, plotted_T6_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 0.01 mu translation, Tz, 5000MU",absolute=False,iternum=9,fig=fig141,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig141.number}.pdf")

fig142=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig142=plot_RMS_bygroup(align_out_T6_5000MU, plotted_T6_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 0.01 mu translation, Rx, 5000MU",absolute=False,iternum=9,fig=fig142,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig142.number}.pdf")

fig143=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig143=plot_RMS_bygroup(align_out_T6_5000MU, plotted_T6_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 0.01 mu translation, Rz, 5000MU",absolute=False,iternum=9,fig=fig143,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig143.number}.pdf")

# 7
fig144=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig144=plot_RMS_bygroup(align_out_T7_5000MU, plotted_T7_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 0.001 mu translation, Tx, 5000MU",absolute=False,iternum=9,fig=fig144,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig144.number}.pdf")

fig145=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig145=plot_RMS_bygroup(align_out_T7_5000MU, plotted_T7_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 0.001 mu translation, Tz, 5000MU",absolute=False,iternum=9,fig=fig145,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig145.number}.pdf")

fig146=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig146=plot_RMS_bygroup(align_out_T7_5000MU, plotted_T7_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 0.001 mu translation, Rx, 5000MU",absolute=False,iternum=9,fig=fig146,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig146.number}.pdf")

fig147=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig147=plot_RMS_bygroup(align_out_T7_5000MU, plotted_T7_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 0.001 mu translation, Rz, 5000MU",absolute=False,iternum=9,fig=fig147,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig147.number}.pdf")


### high momentum tracks
# 0
fig148=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig148=plot_RMS_bygroup(align_out_HT0_5000MU, plotted_HT0_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 200 mu HighETracks, Tx, 5000MU",absolute=False,iternum=9,fig=fig148,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig148.number}.pdf")

fig149=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig149=plot_RMS_bygroup(align_out_HT0_5000MU, plotted_HT0_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 200 mu HighETracks, Tz, 5000MU",absolute=False,iternum=9,fig=fig149,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig149.number}.pdf")

fig150=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig150=plot_RMS_bygroup(align_out_HT0_5000MU, plotted_HT0_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 200 mu HighETracks, Rx, 5000MU",absolute=False,iternum=9,fig=fig150,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig150.number}.pdf")

fig151=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig151=plot_RMS_bygroup(align_out_HT0_5000MU, plotted_HT0_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 200 mu HighETracks, Rz, 5000MU",absolute=False,iternum=9,fig=fig151,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig151.number}.pdf")

# 1
fig152=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig152=plot_RMS_bygroup(align_out_HT1_5000MU, plotted_HT1_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 200 mu HighETracks, Tx, 5000MU",absolute=False,iternum=9,fig=fig152,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig152.number}.pdf")

fig153=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig153=plot_RMS_bygroup(align_out_HT1_5000MU, plotted_HT1_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 200 mu HighETracks, Tz, 5000MU",absolute=False,iternum=9,fig=fig153,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig153.number}.pdf")

fig154=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig154=plot_RMS_bygroup(align_out_HT1_5000MU, plotted_HT1_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 200 mu HighETracks, Rx, 5000MU",absolute=False,iternum=9,fig=fig154,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig154.number}.pdf")

fig155=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig155=plot_RMS_bygroup(align_out_HT1_5000MU, plotted_HT1_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 200 mu HighETracks, Rz, 5000MU",absolute=False,iternum=9,fig=fig155,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig155.number}.pdf")

# 2
fig156=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig156=plot_RMS_bygroup(align_out_HT2_5000MU, plotted_HT2_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 200 mu HighETracks, Tx, 5000MU",absolute=False,iternum=9,fig=fig156,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig156.number}.pdf")

fig157=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig157=plot_RMS_bygroup(align_out_HT2_5000MU, plotted_HT2_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 200 mu HighETracks, Tz, 5000MU",absolute=False,iternum=9,fig=fig157,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig157.number}.pdf")

fig158=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig158=plot_RMS_bygroup(align_out_HT2_5000MU, plotted_HT2_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 200 mu HighETracks, Rx, 5000MU",absolute=False,iternum=9,fig=fig158,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig158.number}.pdf")

fig159=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig159=plot_RMS_bygroup(align_out_HT2_5000MU, plotted_HT2_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 200 mu HighETracks, Rz, 5000MU",absolute=False,iternum=9,fig=fig159,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig159.number}.pdf")

# 3
fig160=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig160=plot_RMS_bygroup(align_out_HT3_5000MU, plotted_HT3_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 200 mu HighETracks, Tx, 5000MU",absolute=False,iternum=9,fig=fig160,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig160.number}.pdf")

fig161=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig161=plot_RMS_bygroup(align_out_HT3_5000MU, plotted_HT3_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 200 mu HighETracks, Tz, 5000MU",absolute=False,iternum=9,fig=fig161,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig161.number}.pdf")

fig162=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig162=plot_RMS_bygroup(align_out_HT3_5000MU, plotted_HT3_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 200 mu HighETracks, Rx, 5000MU",absolute=False,iternum=9,fig=fig162,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig162.number}.pdf")

fig163=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig163=plot_RMS_bygroup(align_out_HT3_5000MU, plotted_HT3_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 200 mu HighETracks, Rz, 5000MU",absolute=False,iternum=9,fig=fig163,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig163.number}.pdf")

# 4
fig164=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig164=plot_RMS_bygroup(align_out_HT4_5000MU, plotted_HT4_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 200 mu HighETracks, Tx, 5000MU",absolute=False,iternum=9,fig=fig164,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig164.number}.pdf")

fig165=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig165=plot_RMS_bygroup(align_out_HT4_5000MU, plotted_HT4_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 200 mu HighETracks, Tz, 5000MU",absolute=False,iternum=9,fig=fig165,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig165.number}.pdf")

fig166=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig166=plot_RMS_bygroup(align_out_HT4_5000MU, plotted_HT4_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 200 mu HighETracks, Rx, 5000MU",absolute=False,iternum=9,fig=fig166,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig166.number}.pdf")

fig167=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig167=plot_RMS_bygroup(align_out_HT4_5000MU, plotted_HT4_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 200 mu HighETracks, Rz, 5000MU",absolute=False,iternum=9,fig=fig167,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig167.number}.pdf")

# 5
fig168=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig168=plot_RMS_bygroup(align_out_HT5_5000MU, plotted_HT5_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 200 mu HighETracks, Tx, 5000MU",absolute=False,iternum=9,fig=fig168,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig168.number}.pdf")

fig169=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig169=plot_RMS_bygroup(align_out_HT5_5000MU, plotted_HT5_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 200 mu HighETracks, Tz, 5000MU",absolute=False,iternum=9,fig=fig169,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig169.number}.pdf")

fig170=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig170=plot_RMS_bygroup(align_out_HT5_5000MU, plotted_HT5_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 200 mu HighETracks, Rx, 5000MU",absolute=False,iternum=9,fig=fig170,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig170.number}.pdf")

fig171=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig171=plot_RMS_bygroup(align_out_HT5_5000MU, plotted_HT5_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 200 mu HighETracks, Rz, 5000MU",absolute=False,iternum=9,fig=fig171,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig171.number}.pdf")

# 6
fig172=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig172=plot_RMS_bygroup(align_out_HT6_5000MU, plotted_HT6_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 200 mu HighETracks, Tx, 5000MU",absolute=False,iternum=9,fig=fig172,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig172.number}.pdf")

fig173=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig173=plot_RMS_bygroup(align_out_HT6_5000MU, plotted_HT6_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 200 mu HighETracks, Tz, 5000MU",absolute=False,iternum=9,fig=fig173,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig173.number}.pdf")

fig174=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig174=plot_RMS_bygroup(align_out_HT6_5000MU, plotted_HT6_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 200 mu HighETracks, Rx, 5000MU",absolute=False,iternum=9,fig=fig174,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig174.number}.pdf")

fig175=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig175=plot_RMS_bygroup(align_out_HT6_5000MU, plotted_HT6_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 200 mu HighETracks, Rz, 5000MU",absolute=False,iternum=9,fig=fig175,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig175.number}.pdf")

# 7
fig176=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig176=plot_RMS_bygroup(align_out_HT7_5000MU, plotted_HT7_5000MU,title_t,filetitle,"Tx",newgrouping,"misalign 200 mu HighETracks, Tx, 5000MU",absolute=False,iternum=9,fig=fig176,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig176.number}.pdf")

fig177=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig177=plot_RMS_bygroup(align_out_HT7_5000MU, plotted_HT7_5000MU,title_t,filetitle,"Tz",newgrouping,"misalign 200 mu HighETracks, Tz, 5000MU",absolute=False,iternum=9,fig=fig177,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig177.number}.pdf")

fig178=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig178=plot_RMS_bygroup(align_out_HT7_5000MU, plotted_HT7_5000MU,title_r,filetitle,"Rx",newgrouping,"misalogn 200 mu HighETracks, Rx, 5000MU",absolute=False,iternum=9,fig=fig178,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig178.number}.pdf")

fig179=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config 5, 2020 data, 5000MU",absolute=False,iternum=9)
fig179=plot_RMS_bygroup(align_out_HT7_5000MU, plotted_HT7_5000MU,title_r,filetitle,"Rz",newgrouping,"misalign 200 mu HighETracks, Rz, 5000MU",absolute=False,iternum=9,fig=fig179,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig179.number}.pdf")

### new combi 1
#fig1=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_t,filetitle,"Tx",newgrouping,"Tx, flo, 5000MU",absolute=False,iternum=9)
#fig1=plot_RMS_bygroup(align_out_conf5_Rz_5000MU, plotted_conf5_Rz_5000MU,title_t,filetitle,"Tx",newgrouping,"Tx, combi5 + backlayer Rz, 5000MU",absolute=False,iternum=9,fig=fig1,color="C1",nudge=20)
#plt.ylabel(r"Mean Tx ($\mu m$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_conf5_Rz_Tx_5000MU_{fig1.number}.pdf")

#fig2=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_t,filetitle,"Tz",newgrouping,"fullL, Tz, flo, 5000MU",absolute=False,iternum=9)
#fig2=plot_RMS_bygroup(align_out_conf5_Rz_5000MU, plotted_conf5_Rz_5000MU,title_t,filetitle,"Tz",newgrouping,"Tz, combi5 + backlayer Rz, 5000MU",absolute=False,iternum=9,fig=fig2,color="C1",nudge=20)
#plt.ylabel(r"Mean Tz ($\mu m$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_conf5_Rz_Tz_5000MU_{fig2.number}.pdf")

#fig3=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_r,filetitle,"Rx",newgrouping,"fullL, Rx, flo, 5000MU",absolute=False,iternum=9)
#fig3=plot_RMS_bygroup(align_out_conf5_Rz_5000MU, plotted_conf5_Rz_5000MU,title_r,filetitle,"Rx",newgrouping,"Rx, combi5 + backlayer Rz, 5000MU",absolute=False,iternum=9,fig=fig3,color="C1",nudge=20)
#plt.ylabel(r"Mean Rx (mrad)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_conf5_Rz_Rx_5000MU_{fig3.number}.pdf")

#fig4=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_r,filetitle,"Rz",newgrouping,"fullL, Rz, flo, 5000MU",absolute=False,iternum=9)
#fig4=plot_RMS_bygroup(align_out_conf5_Rz_5000MU, plotted_conf5_Rz_5000MU,title_r,filetitle,"Rz",newgrouping,"Rz, combi5 + backlayer Rz, 5000MU",absolute=False,iternum=9,fig=fig4,color="C1",nudge=20)
#plt.ylabel(r"Mean Rz (mrad$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_conf5_Rz_Rz_5000MU_{fig4.number}.pdf")

### new combi 2
#fig5=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_t,filetitle,"Tx",newgrouping,"fullL, Tx, flo, 5000MU",absolute=False,iternum=9)
#fig5=plot_RMS_bygroup(align_out_noTx_5000MU, plotted_noTx_5000MU,title_t,filetitle,"Tx",newgrouping,"noTx in backlayer, Tx, 5000MU",absolute=False,iternum=9,fig=fig5,color="C1",nudge=20)
#plt.ylabel(r"Mean Tx ($\mu m$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_noTx_Tx_5000MU_{fig5.number}.pdf")

#fig6=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_t,filetitle,"Tz",newgrouping,"fullL, Tz, flo, 5000MU",absolute=False,iternum=9)
#fig6=plot_RMS_bygroup(align_out_noTx_5000MU, plotted_noTx_5000MU,title_t,filetitle,"Tz",newgrouping,"noTx in backlayer, Tz, 5000MU",absolute=False,iternum=9,fig=fig6,color="C1",nudge=20)
#plt.ylabel(r"Mean Tz ($\mu m$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_noTx_Tz_5000MU_{fig6.number}.pdf")

#fig7=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_r,filetitle,"Rx",newgrouping,"fullL, Rx, flo, 5000MU",absolute=False,iternum=9)
#fig7=plot_RMS_bygroup(align_out_noTx_5000MU, plotted_noTx_5000MU,title_r,filetitle,"Rx",newgrouping,"noTx in backlayer, Rx, 5000MU",absolute=False,iternum=9,fig=fig7,color="C1",nudge=20)
#plt.ylabel(r"Mean Rx (mrad)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_noTx_Rx_5000MU_{fig7.number}.pdf")

#fig8=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_r,filetitle,"Rz",newgrouping,"fullL, Rz, flo, 5000MU",absolute=False,iternum=9)
#fig8=plot_RMS_bygroup(align_out_noTx_5000MU, plotted_noTx_5000MU,title_r,filetitle,"Rz",newgrouping,"noTx in backlayer, Rz, 5000MU",absolute=False,iternum=9,fig=fig8,color="C1",nudge=20)
#plt.ylabel(r"Mean Rz (mrad)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_noTx_Rz_5000MU_{fig8.number}.pdf")

## TxRz on config 5
#fig9=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_t,filetitle,"Tx",newgrouping,"fullL, Tx, flo, 5000MU",absolute=False,iternum=9)
#fig9=plot_RMS_bygroup(align_out_TxRz_5000MU, plotted_TxRz_5000MU,title_t,filetitle,"Tx",newgrouping,"TxRz in backlayer, Tx, 5000MU",absolute=False,iternum=9,fig=fig9,color="C1",nudge=20)
#plt.ylabel(r"Mean Tx ($\mu m$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_TxRz_Tx_5000MU_{fig9.number}.pdf")

#fig10=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_t,filetitle,"Tz",newgrouping,"fullL, Tz, flo, 5000MU",absolute=False,iternum=9)
#fig10=plot_RMS_bygroup(align_out_TxRz_5000MU, plotted_TxRz_5000MU,title_t,filetitle,"Tz",newgrouping,"TxRz in backlayer, Tz, 5000MU",absolute=False,iternum=9,fig=fig10,color="C1",nudge=20)
#plt.ylabel(r"Mean Tz ($\mu m$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_TxRz_Tz_5000MU_{fig10.number}.pdf")

#fig11=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_r,filetitle,"Rx",newgrouping,"fullL, Rx, flo, 5000MU",absolute=False,iternum=9)
#fig11=plot_RMS_bygroup(align_out_TxRz_5000MU, plotted_TxRz_5000MU,title_r,filetitle,"Rx",newgrouping,"TxRz in backlayer, Rx, 5000MU",absolute=False,iternum=9,fig=fig11,color="C1",nudge=20)
#plt.ylabel(r"Mean Rx (mrad)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_TxRz_Rx_5000MU_{fig11.number}.pdf")

#fig12=plot_RMS_bygroup(align_out_Flo12_5000MU, plotted_Flo12_5000MU, title_r,filetitle,"Rz",newgrouping,"fullL, Rz, flo, 5000MU",absolute=False,iternum=9)
#fig12=plot_RMS_bygroup(align_out_TxRz_5000MU, plotted_TxRz_5000MU,title_r,filetitle,"Rz",newgrouping,"TxRz in backlayer, Rz, 5000MU",absolute=False,iternum=9,fig=fig12,color="C1",nudge=20)
#plt.ylabel(r"Mean Rz (mrad)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/new_TxRz_Rz_5000MU_{fig12.number}.pdf")

################# 2019 + 2020 data ##################
# MU 2020
fig13=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config5, 2020 data, 5000MU",absolute=False,iternum=9)
fig13=plot_RMS_bygroup(align_out_2020_5000MU, plotted_2020_5000MU,title_t,filetitle,"Tx",newgrouping,"2020 data, config5 with Rz, Tx, 5000MU",absolute=False,iternum=9,fig=fig13,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2020data_Tx_5000MU_{fig13.number}.pdf")

fig14=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config5, 2020 data, 5000MU",absolute=False,iternum=9)
fig14=plot_RMS_bygroup(align_out_2020_5000MU, plotted_2020_5000MU,title_t,filetitle,"Tz",newgrouping,"2020 data, config5 with Rz, Tz, 5000MU",absolute=False,iternum=9,fig=fig14,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2020data_Tz_5000MU_{fig14.number}.pdf")

fig15=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config5, 2020 data, 5000MU",absolute=False,iternum=9)
fig15=plot_RMS_bygroup(align_out_2020_5000MU, plotted_2020_5000MU,title_r,filetitle,"Rx",newgrouping,"2020 data, config5 with Rz, Rx, 5000MU",absolute=False,iternum=9,fig=fig15,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2020data_Rx_5000MU_{fig15.number}.pdf")

fig16=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config5, 2020 data, 5000MU",absolute=False,iternum=9)
fig16=plot_RMS_bygroup(align_out_2020_5000MU, plotted_2020_5000MU,title_r,filetitle,"Rz",newgrouping,"2020 data, config5 with Rz, Rz, 5000MU",absolute=False,iternum=9,fig=fig16,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2020data_Rz_5000MU_{fig16.number}.pdf")

# MD 2019
fig17=plot_RMS_bygroup(align_out_config5_5000MD_2019, plotted_config5_5000MD_2019, title_t,filetitle,"Tx",newgrouping,"config5, 2019 data, 5000MD",absolute=False,iternum=9)
fig17=plot_RMS_bygroup(align_out_2019_5000MD, plotted_2019_5000MD,title_t,filetitle,"Tx",newgrouping,"2019 data, config5, Tx, 5000MD",absolute=False,iternum=9,fig=fig17,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2019data_Tx_5000MD_{fig17.number}.pdf")

fig18=plot_RMS_bygroup(align_out_config5_5000MD_2019, plotted_config5_5000MD_2019, title_t,filetitle,"Tz",newgrouping,"config5, 2019 data, 5000MD",absolute=False,iternum=9)
fig18=plot_RMS_bygroup(align_out_2019_5000MD, plotted_2019_5000MD,title_t,filetitle,"Tz",newgrouping,"2019 data, config5, TTz, 5000MD",absolute=False,iternum=9,fig=fig18,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2019data_Tz_5000MD_{fig18.number}.pdf")

fig19=plot_RMS_bygroup(align_out_config5_5000MD_2019, plotted_config5_5000MD_2019, title_r,filetitle,"Rx",newgrouping,"config5, 2019 data, 5000MD",absolute=False,iternum=9)
fig19=plot_RMS_bygroup(align_out_2019_5000MD, plotted_2019_5000MD,title_r,filetitle,"Rx",newgrouping,"2019 data, config5, TRx, 5000MD",absolute=False,iternum=9,fig=fig19,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2019data_Rx_5000MD_{fig19.number}.pdf")

fig20=plot_RMS_bygroup(align_out_config5_5000MD_2019, plotted_config5_5000MD_2019, title_r,filetitle,"Rz",newgrouping,"config5, 2019 data, 5000MD",absolute=False,iternum=9)
fig20=plot_RMS_bygroup(align_out_2019_5000MD, plotted_2019_5000MD,title_r,filetitle,"Rz",newgrouping,"2019 data, config5, TRz, 5000MD",absolute=False,iternum=9,fig=fig20,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2019data_Rz_5000MD_{fig20.number}.pdf")

# MD 2020
fig21=plot_RMS_bygroup(align_out_config5_5000MD_2020, plotted_config5_5000MD_2020, title_t,filetitle,"Tx",newgrouping,"config5, 2020 data, 5000MD",absolute=False,iternum=9)
fig21=plot_RMS_bygroup(align_out_2020_5000MD, plotted_2020_5000MD,title_t,filetitle,"Tx",newgrouping,"2020 data, config5, Tx, 5000MD",absolute=False,iternum=9,fig=fig21,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2020data_Tx_5000MD_{fig21.number}.pdf")

fig22=plot_RMS_bygroup(align_out_config5_5000MD_2020, plotted_config5_5000MD_2020, title_t,filetitle,"Tz",newgrouping,"config5, 2020 data, 5000MD",absolute=False,iternum=9)
fig22=plot_RMS_bygroup(align_out_2020_5000MD, plotted_2020_5000MD,title_t,filetitle,"Tz",newgrouping,"2020 data, config5, Tz, 5000MD",absolute=False,iternum=9,fig=fig22,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2020data_Tz_5000MD_{fig22.number}.pdf")

fig23=plot_RMS_bygroup(align_out_config5_5000MD_2020, plotted_config5_5000MD_2020, title_r,filetitle,"Rx",newgrouping,"config5, 2020 data, 5000MD",absolute=False,iternum=9)
fig23=plot_RMS_bygroup(align_out_2020_5000MD, plotted_2020_5000MD,title_r,filetitle,"Rx",newgrouping,"2020 data, config5, Rx, 5000MD",absolute=False,iternum=9,fig=fig23,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2020data_Rx_5000MD_{fig23.number}.pdf")

fig24=plot_RMS_bygroup(align_out_config5_5000MD_2020, plotted_config5_5000MD_2020, title_r,filetitle,"Rz",newgrouping,"config5, 2020 data, 5000MD",absolute=False,iternum=9)
fig24=plot_RMS_bygroup(align_out_2020_5000MD, plotted_2020_5000MD,title_r,filetitle,"Rz",newgrouping,"2020 data, config5, Rz, 5000MD",absolute=False,iternum=9,fig=fig24,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2020data_Rz_5000MD_{fig24.number}.pdf")

# MU 2019
fig25=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tx",newgrouping,"config5, 2019 data, 5000MU",absolute=False,iternum=9)
fig25=plot_RMS_bygroup(align_out_2019_5000MU, plotted_2019_5000MU,title_t,filetitle,"Tx",newgrouping,"2019 data, config5 with Rz, Tx, 5000MU",absolute=False,iternum=9,fig=fig25,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2019data_Tx_5000MU_{fig25.number}.pdf")

fig26=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_t,filetitle,"Tz",newgrouping,"config5, 2019 data, 5000MU",absolute=False,iternum=9)
fig26=plot_RMS_bygroup(align_out_2019_5000MU, plotted_2019_5000MU,title_t,filetitle,"Tz",newgrouping,"2019 data, config5 with Rz, Tz, 5000MU",absolute=False,iternum=9,fig=fig26,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2019data_Tz_5000MU_{fig26.number}.pdf")

fig27=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rx",newgrouping,"config5, 2019 data, 5000MU",absolute=False,iternum=9)
fig27=plot_RMS_bygroup(align_out_2019_5000MU, plotted_2019_5000MU,title_r,filetitle,"Rx",newgrouping,"2019 data, config5 with Rz, Rx, 5000MU",absolute=False,iternum=9,fig=fig27,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2019data_Rx_5000MU_{fig27.number}.pdf")

fig28=plot_RMS_bygroup(align_out_config5_5000MU, plotted_config5_5000MU, title_r,filetitle,"Rz",newgrouping,"config5, 2019 data, 5000MU",absolute=False,iternum=9)
fig28=plot_RMS_bygroup(align_out_2019_5000MU, plotted_2019_5000MU,title_r,filetitle,"Rz",newgrouping,"2019 data, config5 with Rz, Rz, 5000MU",absolute=False,iternum=9,fig=fig28,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/2019data_Rz_5000MU_{fig28.number}.pdf")

## plotte MU gegen MD 2020
fig29=plot_RMS_bygroup(align_out_2020_5000MU, plotted_2020_5000MU, title_t,filetitle,"Tx",newgrouping,"config5 with Rz, 2020 data, 5000MU",absolute=False,iternum=9)
fig29=plot_RMS_bygroup(align_out_2020_5000MD, plotted_2020_5000MD,title_t,filetitle,"Tx",newgrouping,"2020 data, config5 with Rz, Tx, 5000MD",absolute=False,iternum=9,fig=fig29,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/MUvsMD_Tx_5000MU_{fig29.number}.pdf")

fig30=plot_RMS_bygroup(align_out_2020_5000MU, plotted_2020_5000MU, title_t,filetitle,"Tz",newgrouping,"config5 with Rz, 2020 data, 5000MU",absolute=False,iternum=9)
fig30=plot_RMS_bygroup(align_out_2020_5000MD, plotted_2020_5000MD,title_t,filetitle,"Tz",newgrouping,"2020 data, config5 with Rz, Tz, 5000MD",absolute=False,iternum=9,fig=fig30,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/MUvsMD_Tz_5000MU_{fig30.number}.pdf")

fig31=plot_RMS_bygroup(align_out_2020_5000MU, plotted_2020_5000MU, title_r,filetitle,"Rx",newgrouping,"config5 with Rz, 2020 data, 5000MU",absolute=False,iternum=9)
fig31=plot_RMS_bygroup(align_out_2020_5000MD, plotted_2020_5000MD,title_r,filetitle,"Rx",newgrouping,"2020 data, config5 with Rz, Rx, 5000MD",absolute=False,iternum=9,fig=fig31,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/MUvsMD_Rx_5000MU_{fig31.number}.pdf")

fig32=plot_RMS_bygroup(align_out_2020_5000MU, plotted_2020_5000MU, title_r,filetitle,"Rz",newgrouping,"config5 with Rz, 2020 data, 5000MU",absolute=False,iternum=9)
fig32=plot_RMS_bygroup(align_out_2020_5000MD, plotted_2020_5000MD,title_r,filetitle,"Rz",newgrouping,"2020 data, config5 with Rz, Tz, 5000MD",absolute=False,iternum=9,fig=fig32,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/MUvsMD_Rz_5000MU_{fig32.number}.pdf")

########## and for 2019
fig33=plot_RMS_bygroup(align_out_2019_5000MU, plotted_2019_5000MU, title_t,filetitle,"Tx",newgrouping,"config5 with Rz, 2019 data, 5000MU",absolute=False,iternum=9)
fig33=plot_RMS_bygroup(align_out_2019_5000MD, plotted_2019_5000MD,title_t,filetitle,"Tx",newgrouping,"2019 data, config5 with Rz, Tx, 5000MD",absolute=False,iternum=9,fig=fig33,color="C1",nudge=20)
plt.ylabel(r"Mean Tx ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/MUvsMD_Tx_5000MU_{fig33.number}.pdf")

fig34=plot_RMS_bygroup(align_out_2019_5000MU, plotted_2019_5000MU, title_t,filetitle,"Tz",newgrouping,"config5 with Rz, 2019 data, 5000MU",absolute=False,iternum=9)
fig34=plot_RMS_bygroup(align_out_2019_5000MD, plotted_2019_5000MD,title_t,filetitle,"Tz",newgrouping,"2019 data, config5 with Rz, Tz, 5000MD",absolute=False,iternum=9,fig=fig34,color="C1",nudge=20)
plt.ylabel(r"Mean Tz ($\mu m$)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/MUvsMD_Tz_5000MU_{fig34.number}.pdf")

fig35=plot_RMS_bygroup(align_out_2019_5000MU, plotted_2019_5000MU, title_r,filetitle,"Rx",newgrouping,"config5 with Rz, 2019 data, 5000MU",absolute=False,iternum=9)
fig35=plot_RMS_bygroup(align_out_2019_5000MD, plotted_2019_5000MD,title_r,filetitle,"Rx",newgrouping,"2019 data, config5 with Rz, Rx, 5000MD",absolute=False,iternum=9,fig=fig35,color="C1",nudge=20)
plt.ylabel(r"Mean Rx (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/MUvsMD_Rx_5000MU_{fig35.number}.pdf")

fig36=plot_RMS_bygroup(align_out_2019_5000MU, plotted_2019_5000MU, title_r,filetitle,"Rz",newgrouping,"config5 with Rz, 2019 data, 5000MU",absolute=False,iternum=9)
fig36=plot_RMS_bygroup(align_out_2019_5000MD, plotted_2019_5000MD,title_r,filetitle,"Rz",newgrouping,"2019 data, config5 with Rz, Tz, 5000MD",absolute=False,iternum=9,fig=fig36,color="C1",nudge=20)
plt.ylabel(r"Mean Rz (mrad)")
plt.xlabel("Group position in z (m)")
plt.legend(loc='upper center')
plt.savefig(f"{outdir}/MUvsMD_Rz_5000MU_{fig36.number}.pdf")

##########################

plt.show()
