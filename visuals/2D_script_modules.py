import argparse

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import copy

import json
import numpy as np

## Collecting information from json file
parser = argparse.ArgumentParser()
parser.add_argument("filename_in")
parser.add_argument("foldername_out")
parser.add_argument("-framelayers",action="store_true")
args=parser.parse_args()
thisfile=args.filename_in
outdir=args.foldername_out

#thisfile="20210128_Modules_Tall/parsed_log.json"
#outdir="20210128_Modules_Tall/"

with open(thisfile) as f:
    align_output=json.load(f)

convergences=align_output.pop("converged")
listConvergences=[index for (index,value) in enumerate(convergences) if value==1]
if listConvergences:
    firstConverged=listConvergences[0]
else:
    firstConverged=-1

#fix floats

align_output.pop("total_chi2_vals", None)
align_output.pop("total_chi2_dofs", None)

labels=["Tx","Tz"]
for alignable in align_output.keys():
    for label in labels:
        align_output[alignable][label]=[float(ele.strip(',')) for ele in align_output[alignable][label]]

# Add cone degeneracy test variable
#newlabel="ConeDegeneracy"
#for alignable in align_output.keys():
#    Txlist=align_output[alignable]["Tx"]
#    Tzlist=align_output[alignable]["Tz"]
#    print("---------", alignable, "---------")
#    if Txlist[0] < 1e-8 or Tzlist[0]< 1e-8:
#        continue
#    print(Txlist)
#    print(Tzlist)
#    newlist=[Tx/Tz for Tx,Tz in zip(Txlist,Tzlist)]
#    align_output[alignable][newlabel]=newlist
#labels.append(newlabel)
        
#for non-pandas build: need to remap format of dictionary
alignables=align_output.keys()

if args.framelayers:
    typelabel="Framelayers"
    alignables=[alignable for alignable in alignables if "Side" in alignable]
    
    layertest={}
    for label in labels:
        layertest[label]={}
        for alignable in alignables:
            intro=alignable.split("Side")[0][0:-1]
            quarterindicator=alignable.split("Side")[0][-1]
            
            if "T3" in intro:
                modules=[0,1,2,3,4,5]
            else:
                modules=[0,1,2,3,4]
            for module in modules:
                newalignablefirst=f"{intro}Quarter{'0' if 'C' in quarterindicator else '1'}Module{module}"
                newalignablesecond=f"{intro}Quarter{'2' if 'C' in quarterindicator else '3'}Module{module}"
                
                layertest[label][newalignablefirst]=1000*align_output[alignable][label][-1]
                layertest[label][newalignablesecond]=1000*align_output[alignable][label][-1]
else:
    typelabel="Modules"
    # remove stations, layers from consideration
    alignables=[alignable for alignable in alignables if "Module" in alignable]
    layertest={}
    for label in labels:
        layertest[label]={}
        for alignable in alignables:
            intro=alignable.split("Quarter")[0]
            quarters=alignable.split("|")
            module=alignable.split("Module")

            #split quarters for drawing
            newalignablefirst=f"{intro}Quarter{quarters[0][-1]}Module{module[1][0]}"
            newalignablesecond=f"{intro}Quarter{quarters[1][0]}Module{module[1][0]}"

            # Plot position changes in micrometres
            if "T" in label:
                layertest[label][newalignablefirst]=1000*align_output[alignable][label][-1]
                layertest[label][newalignablesecond]=1000*align_output[alignable][label][-1]
            else:
                layertest[label][newalignablefirst]=align_output[alignable][label][-1]
                layertest[label][newalignablesecond]=align_output[alignable][label][-1]
        
## Setup classes for SciFi modules/mats plots
class scifi_all:
    x_bins, y_bins = 12, 2
    x_start, x_end = -30, 30
    y_start, y_end = -20, 20
    extent = [x_start, x_end, y_start, y_end]

    def __init__(self):
        self.cmax = self.cmin = None
        self.colormap = plt.get_cmap('inferno')
        self.textcolor = 'teal'
        self.p_modules = []
        self.p_mats = []
        self.map_alignable_to_canvas = {}
        #self.map_alignable_to_xzposition = {}

    def reset_modules(self,station=1):
        self.p_modules.clear()
        self.p_mats.clear()
        # Reminder of quarter arrangement
        #  3  |   2
        # ------------
        #  1  |   0

        if station==3:
            modules = [0, 1, 2, 3, 4, 5]
            modulenames = ["Module0","Module1","Module2","Module3","Module4","Module5"]
        else:
            modules = [0,1,2,3,4]
            modulenames= ["Module0","Module1","Module2","Module3","Module4"]
        
        for q in [0, 1, 2, 3]:  # quarters
            for m in modules:  # modules
                self.p_modules.append(patches.Rectangle(
                    (
                        (-5*(m+1) if q in [
                         0, 3] else (5)*m),
                        (0. if q in [0, 1] else scifi_all.y_start)
                    ), 5, scifi_all.y_end, linewidth=3, edgecolor='grey', facecolor='none'))

        self.map_alignable_to_canvas.clear()
        for q, qname in enumerate(["Quarter0", "Quarter1", "Quarter2", "Quarter3"]):  # quarters
            # modules
            for m, mname in enumerate(modulenames):
                self.map_alignable_to_canvas[(
                    0 if q in [2, 3] else 1,
                    5 - m if q in [1, 3] else m + 6
                )] = qname + mname

    def draw_modules(self, ax):
        for m in self.p_modules:
            new_m = copy.copy(m)
            ax.add_patch(new_m)

# modded to accept python dictionary
    def plot_data_modules_singlestation(self, data, title, filename, xlabel='x [a.u.]', ylabel='y [a.u.]', zlabel='',
                                        layerinfo=["FT/T1LayerX1","FT/T1LayerU","FT/T1LayerV","FT/T1LayerX2"], save=True):
        fig, axes = plt.subplots(nrows=1,ncols=len(layerinfo), figsize=(6*len(layerinfo)+1, 4))
#        ax.autoscale()

        self.reset_modules()
        for layer,ax in zip(layerinfo,axes.flat):
            matrix_data = np.zeros((2, 12))
            jump_x = (scifi_all.x_end - scifi_all.x_start) / \
                (2.*scifi_all.x_bins)
            jump_y = (scifi_all.y_end - scifi_all.y_start) / \
                (2.*scifi_all.y_bins)
            x_positions = np.linspace(
                start=scifi_all.x_start, stop=scifi_all.x_end, num=scifi_all.x_bins, endpoint=False)
            y_positions = np.linspace(
                start=scifi_all.y_start, stop=scifi_all.y_end, num=scifi_all.y_bins, endpoint=False)
            for y_idx, yv in enumerate(y_positions):
                for x_idx, xv in enumerate(x_positions):
                    try:
                        name_alignables = self.map_alignable_to_canvas[(y_idx, x_idx)]
                    except KeyError:
                        print(f"no module at {y_idx} {x_idx}")
                        label=None
                        matrix_data[y_idx,x_idx]=0
                    else:
                        label = data[f'{layer}{name_alignables}']
                    matrix_data[y_idx, x_idx] = label
                    
                    if x_idx>8:
                        print(f'idx:{x_idx}, label:{label}, name:{name_alignables}')

                    text_x = xv + jump_x
                    text_y = yv + jump_y
                    ax.text(text_x, text_y, "%3.2f" % label if label is not None else "", ha='center', va='center',
                            fontsize=12, color=self.textcolor, fontweight='bold', rotation=75)
            im = ax.imshow(matrix_data, cmap=self.colormap, origin='upper',  # viridis
                           extent=scifi_all.extent,
                           vmin=self.cmin, vmax=self.cmax,
                           interpolation=None)
            ax.set_title(f'{title},{layer}' , fontdict={"fontsize": 12})
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            self.draw_modules(ax)
            #ax.set_axis_off()
        fig.subplots_adjust(right=0.8, bottom=0.04)
        cb = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.15, pad=0.02)
        cb.set_label(zlabel)
        
        if save:
            fig.savefig(filename)
            ax.cla()
            fig.clf()
            del fig
            del ax
        else:
            return fig

    def plot_data_modules_multistation(self, data, title, filename, xlabel='x [a.u.]', ylabel='y [a.u.]', zlabel='',
                                        stations=["FT/T1","FT/T2","FT/T3"],layers=["LayerX1","LayerU","LayerV","LayerX2"], save=True):
        fig, axes = plt.subplots(nrows=len(stations),ncols=len(layers), figsize=(5*len(layers)+1, 3*len(stations)))
#        ax.autoscale()
        fig.tight_layout()

        #for layer,ax in zip(layers,axes.flat):
        for row in range(0,len(stations)):
            station=stations[row]
            if "3" in station:
                self.reset_modules(station=3)
            else:
                self.reset_modules()
            for col in range(0,len(layers)):
                layer=layers[col]
                ax=axes[row,col]

                matrix_data = np.zeros((2, 12))
                jump_x = (scifi_all.x_end - scifi_all.x_start) / \
                    (2.*scifi_all.x_bins)
                jump_y = (scifi_all.y_end - scifi_all.y_start) / \
                    (2.*scifi_all.y_bins)
                x_positions = np.linspace(
                    start=scifi_all.x_start, stop=scifi_all.x_end, num=scifi_all.x_bins, endpoint=False)
                y_positions = np.linspace(
                    start=scifi_all.y_start, stop=scifi_all.y_end, num=scifi_all.y_bins, endpoint=False)
                for y_idx, yv in enumerate(y_positions):
                    for x_idx, xv in enumerate(x_positions):
                        try:
                            name_alignables = self.map_alignable_to_canvas[(y_idx, x_idx)]
                        except KeyError:
                            #print(f"no module at {y_idx} {x_idx}")
                            label=None
                            matrix_data[y_idx,x_idx]=0
                        else:
                            try:
                                label = data[f'{station}{layer}{name_alignables}']
                                matrix_data[y_idx, x_idx] = label
                            except KeyError:
                                matrix_data[y_idx, x_idx]=0

                        text_x = xv + jump_x
                        text_y = yv + jump_y
                        ax.text(text_x, text_y, "%3.2f" % label if label is not None else "", ha='center', va='center',
                                fontsize=12, color=self.textcolor, fontweight='bold', rotation=75)
                im = ax.imshow(matrix_data, cmap=self.colormap, origin='upper',  # viridis
                               extent=scifi_all.extent,
                               vmin=self.cmin, vmax=self.cmax,
                               interpolation=None)
                ax.set_title(f'{title},{station}{layer}' , fontdict={"fontsize": 12})
                #ax.set_xlabel(xlabel)
                #ax.set_ylabel(ylabel)
                self.draw_modules(ax)
                #ax.set_axis_off()
                xax = ax.get_xaxis() 
                xax = xax.set_visible(False) 
                yax = ax.get_yaxis() 
                yax = yax.set_visible(False) 
            
        fig.subplots_adjust(right=0.8, bottom=0.04)
        cb = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.15, pad=0.02)
        cb.set_label(zlabel)
        
        if save:
            fig.savefig(filename)
            ax.cla()
            fig.clf()
            del fig
            del ax
        else:
            return fig

### Plot the results from the file
print("###### plot generation ######")
print("..... wait .....")
SF = scifi_all()

values=layertest["Tx"].values()
values=list(values)
mods_max = max( abs(max(values)), abs(min(values)) )
SF.cmax = mods_max
SF.cmin = -mods_max

SF.textcolor = 'black'
SF.colormap = plt.get_cmap('BrBG_r')
result=SF.plot_data_modules_multistation(layertest["Tx"], "Tx",
                     f'{outdir}/{typelabel}-Tx-grid.pdf', zlabel=r'Tx $[\mu m]$',save=True)

values=layertest["Tz"].values()
values=list(values)
mods_max = max( abs(max(values)), abs(min(values)) )
SF.cmax = mods_max
SF.cmin = -mods_max
result=SF.plot_data_modules_multistation(layertest["Tz"], "Tz",
                     f'{outdir}/{typelabel}-Tz-grid.pdf', zlabel=r'Tz $[\mu m]$',save=True)

#mods_max = 3
#SF.cmax = mods_max
#SF.cmin = -mods_max
#result=SF.plot_data_modules_multistation(layertest["ConeDegeneracy"], "ConeDegeneracy (Tx/Tz)",
#                     f'{outdir}/{typelabel}-Cone-grid.pdf', zlabel=r'Cone (Tx/Tz)$',save=True)

print("....plots saved!....")
plt.show()

