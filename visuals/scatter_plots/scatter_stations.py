import argparse

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

#for recolor
import matplotlib.colors as mplcolor
import colorsys as colorsys

import json
import numpy as np


## Collecting information from json file
parser = argparse.ArgumentParser()
parser.add_argument("filename_in")
parser.add_argument("foldername_out")
parser.add_argument("-rotations",action="store_true")
parser.add_argument("-layers",action="store_true")
parser.add_argument("-modules",action="store_true")
parser.add_argument("-mats",action="store_true")
parser.add_argument("-stations", action="store_true")
args=parser.parse_args()
thisfile=args.filename_in
outdir=args.foldername_out

with open(thisfile) as f:
    align_output=json.load(f)

# Configure convergence record)
print(align_output['converged'])
convergences=align_output.pop("converged")
#print(convergences)
firstConverged=next(index for (index,value) in enumerate(convergences) if value==1)
iterations=len(convergences)
#print(firstConverged)

# Configure inputs
if args.rotations:
    labels=["Tx","Ty","Tz","Rx","Ry","Rz"]
else:
    labels=["Tx","Ty","Tz"]
positions=["x_global","y_global","z_global"]
expectations=[0,0,0,0,0,0]

#fix floats
align_output.pop('total_chi2_vals', None)
align_output.pop('total_chi2_dofs', None)

for alignable in align_output.keys():
    for label in labels+positions:
        if label in align_output[alignable].keys():
            align_output[alignable][label]=[float(ele.strip(',')) for ele in align_output[alignable][label]]
        else:
            continue

# Add cone degeneracy test variable (put aay for now since it wont work just yet)
newlabel="ConeDegeneracy"
for alignable in align_output.keys():
    Txlist=align_output[alignable]["Tx"]
    Tzlist=align_output[alignable]["Tz"]
    if 0 in Txlist or 0 in Tzlist:
        print('zero in here')
        continue
    else:
        newlist=[Tx/Tz for Tx,Tz in zip(Txlist,Tzlist)]
        align_output[alignable][newlabel]=newlist
#labels.append(newlabel)
            
# Configure plotting
def scale_lightness(rgb, scale_l):
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

def export_legend(legend, filename="legend.png", expand=[-5,-5,5,5]):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)
    
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
    if "Module0" in inputstring:
        return 0
    elif "Module1" in inputstring:
        return 1
    elif "Module2" in inputstring:
        return 2
    elif "Module3" in inputstring:
        return 3
    elif "Module4" in inputstring:
        return 4
    elif "Module5" in inputstring:
        return 5
    else:
        return -1

def plot_scatter(align_output,plotted_alignables,title,filetitle,colourmap,colourname,labels=labels,positions=positions):
        print("Creating plots per iteration...")
        print('labels: ', labels)
        for label in labels:
            print('current label: ', label)
            plt.figure()
            plt.axhline(0,color='r',label="Expected central value") #from expectations
            plt.title(title)
            plt.xlabel("Iteration count")
            plt.ylabel(f'{label} (mm)' if "T" in label else f'{label} (mrad)' if "R" in label else f'{label}')
            print("--------- 1 -----------")
            for ii,alignable in enumerate(plotted_alignables):
                plt.scatter(np.linspace(0,iterations-1,iterations),
                            align_output[alignable][label],color=colourmap[ii],linestyle="-",marker=".",label=alignable)
                plt.scatter(iterations-1,align_output[alignable][label][-1],color=colourmap[ii],marker="*",label="_nolabel")
            plt.axvline(firstConverged,linestyle=":")
            lgd=plt.legend(bbox_to_anchor=(1.05, 1))
            plt.savefig(f'{outdir}/{filetitle}-{label}-{colourname}.pdf',bbox_inches='tight')
            export_legend(lgd,filename=f"{outdir}/legend-{colourname}.pdf")

        print("Creating plots per 2D combination...")
        for jj,label in enumerate(labels):
            for label2 in labels[jj+1:]:
                plt.figure()
                plt.scatter(0,0,color='r',label="Expected central value")
                plt.title(title)
                plt.xlabel(f'{label} (mm)' if "T" in label else f'{label} (mrad)' if "R" in label else f'{label}')
                plt.ylabel(f'{label2} (mm)' if "T" in label2 else f'{label2} (mrad)' if "R" in label2 else f'{label2}')
                for ii,alignable in enumerate(plotted_alignables):
                    plt.scatter(align_output[alignable][label][-1],align_output[alignable][label2][-1],color=colourmap[ii],marker="*",label="_nolabel")
                plt.savefig(f'{outdir}/{filetitle}-{label}{label2}-{colourname}.pdf',bbox_inches='tight')

            for position in positions:
                plt.figure()
                plt.title(title)
                plt.xlabel(f'{position} (mm)')
                plt.ylabel(f'{label} (mm)' if "T" in label else f'{label} (mrad)' if "R" in label else f'{label}')
                for ii,alignable in enumerate(plotted_alignables):
                    plt.scatter(align_output[alignable][position][-1],
                                align_output[alignable][label][-1],
                                color=colourmap[ii],marker="*",label="_nolabel")
                plt.savefig(f'{outdir}/{filetitle}-{position}{label}-{colourname}.pdf',bbox_inches='tight')
        plt.show()
## Begin plotting sequences:
if args.mats:
    # plots for mats
    print("Mats code not ready!")
if args.modules:
    # Set up list to plot
    plotted_alignables=align_output.keys()
    print(plotted_alignables)
    plotted_alignables = [alignable for alignable in plotted_alignables if "Mat" not in alignable]
    plotted_alignables = [alignable for alignable in plotted_alignables if "Module" in alignable]
    
    # Select colourmap   
    # Default: "technicolour" map
    halfstations_colourmap=list(map(halfstation_to_int,plotted_alignables))
    modules_colourmap=list(map(module_to_int,plotted_alignables))
    colormax=max(halfstations_colourmap)
    modulehalfstations_colormap=[None]*len(modules_colourmap)
    for ii,(module,halfstation) in enumerate(zip(modules_colourmap,halfstations_colourmap)):
        modulehalfstations_colormap[ii]=scale_lightness(mplcolor.to_rgb(f'C{module}'),0.3*halfstation/colormax+(1-0.3))
    #maps scaling into 0.5 to 1.5 range, 1 is no recolor)

    colourmap=modulehalfstations_colormap
    colourname="technicolour"
    
    title="alignment of modules"
    filetitle="Modules_scatter"
    
    plot_scatter(align_output,plotted_alignables,title,filetitle,colourmap,colourname)
    
if args.layers:
    # Set up list to plot
    plotted_alignables=align_output.keys()
    plotted_alignables = [alignable for alignable in plotted_alignables if "Mat" not in alignable]
    plotted_alignables = [alignable for alignable in plotted_alignables if "Module" not in alignable]
    plotted_alignables = [alignable for alignable in plotted_alignables if "Layer" in alignable]
    
    # Select colourmap
    stations_colourmap=[int(alignable[4])-1 for alignable in plotted_alignables] #colour from station number
    colourmap=[f'C{ii}' for ii in stations_colourmap]
    colourname="stationcolour"
    
    title="Alignment of Layers"
    filetitle="Layers_scatter"
    
    plot_scatter(align_output,plotted_alignables,title,filetitle,colourmap,colourname)
    
if args.stations:
    # Set up list to plot
    plotted_alignables=align_output.keys()
    plotted_alignables = [alignable for alignable in plotted_alignables if "Mat" not in alignable]
    plotted_alignables = [alignable for alignable in plotted_alignables if "Module" not in alignable]
    plotted_alignables = [alignable for alignable in plotted_alignables if "Layer" not in alignable]
    # Select colourmap
    stations_colourmap=[int(alignable[4])-1 for alignable in plotted_alignables] #colour from station number
    colourmap=[f'C{ii}' for ii in stations_colourmap]
    colourname="stationcolour"
    
    title="Alignment of Stations"
    filetitle="Stations_scatter"
    plot_scatter(align_output,plotted_alignables,title,filetitle,colourmap,colourname)


## data ##
# florian changes, without constraints, e = 1000
#thisfile = "../output_alignment/20212605/MU/n10_e1000/updateFlorian/AlignmentResults/parsedlog.json"
#align_out_Flo10_1000MU = open_alignment(thisfile)
#plotted_Flo10_1000MU = choose_modules(align_out_Flo10_1000MU)

# florian changes without constraints, e = 3000
#thisfile = "../output_alignment/20212605/MU/n10_e3000/updateFlorian/AlignmentResults/parsedlog.json"
#align_out_Flo10_3000MU = open_alignment(thisfile)
#plotted_Flo10_3000MU = choose_modules(align_out_Flo10_3000MU)

#thisfile = "../output_alignment/20211406/MU/n10_e1000/T3_allT_freeRx/AlignmentResults/parsedlog.json"
#align_out_T3FreeRx_1000MU = open_alignment(thisfile)
#plotted_T3FreeRx_1000MU = choose_modules(align_out_T3FreeRx_1000MU)
## data end ##

## start plotting

#fig1=plot_RMS_bygroup(align_out_Flo10_1000MU, plotted_Flo10_1000MU, title,filetitle,"Tx",newgrouping,"fullL, Tx, flo, 1000MU",absolute=False,iternum=9)
#fig1=plot_RMS_bygroup(align_out_T3FreeRx_1000MU, plotted_T3FreeRx_1000MU,title,filetitle,"Tx",newgrouping,"fullL, Tx, T3 Rx Free, 1000MU",absolute=False,iternum=9,fig=fig1,color="C1",nudge=20)
#plt.ylabel(r"Mean Tx ($\mu m$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/T3FreeRx_Tx_1000Mu_{fig1.number}.pdf")

#fig2=plot_RMS_bygroup(align_out_Flo10_1000MU, plotted_Flo10_1000MU, title,filetitle,"Tz",newgrouping,"fullL, Tz, flo, 1000MU",absolute=False,iternum=9)
#fig2=plot_RMS_bygroup(align_out_T3FreeRx_1000MU, plotted_T3FreeRx_1000MU,title,filetitle,"Tz",newgrouping,"fullL, Tz, T3 Rx Free, 1000MU",absolute=False,iternum=9,fig=fig2,color="C1",nudge=20)
#plt.ylabel(r"Mean Tz ($\mu m$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/T3FreeRx_Tz_1000Mu_{fig2.number}.pdf")

#fig3=plot_RMS_bygroup(align_out_Flo10_1000MU, plotted_Flo10_1000MU, title,filetitle,"Rx",newgrouping,"fullL, Rx, flo, 1000MU",absolute=False,iternum=9)
#fig3=plot_RMS_bygroup(align_out_T3FreeRx_1000MU, plotted_T3FreeRx_1000MU,title,filetitle,"Rx",newgrouping,"fullL, Rx, T3 Rx Free, 1000MU",absolute=False,iternum=9,fig=fig3,color="C1",nudge=20)
#plt.ylabel(r"Mean Rx ($mrad$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/T3FreeRx_Rx_1000Mu_{fig3.number}.pdf")

#fig4=plot_RMS_bygroup(align_out_Flo10_1000MU, plotted_Flo10_1000MU, title,filetitle,"Rz",newgrouping,"fullL, Rz, flo, 1000MU",absolute=False,iternum=9)
#fig4=plot_RMS_bygroup(align_out_T3FreeRx_1000MU, plotted_T3FreeRx_1000MU,title,filetitle,"Rz",newgrouping,"fullL, Rz, T3 Rx Free, 1000MU",absolute=False,iternum=9,fig=fig4,color="C1",nudge=20)
#plt.ylabel(r"Mean Rz ($mrad$)")
#plt.xlabel("Group position in z (m)")
#plt.legend(loc='upper center')
#plt.savefig(f"{outdir}/T3FreeRx_Rz_1000Mu_{fig4.number}.pdf")

## end plotting
