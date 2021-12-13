import matplotlib.pyplot as plt
import matplotlib.patches as patches
import copy

#for recolor
import matplotlib.colors as mplcolor
import colorsys as colorsys

import json
import numpy as np

thisfile="../output/alignment_runs/done_runs/HE_runs/9/AlignmentResults/parsedlog.json"
outdir="../output/alignment_runs/done_runs/HE_runs/9/colorplots/"
#thisfile="20210412_Modules_allT-10000MU/AlignmentResults/parsed_log.json"
#outdir="20210412_Modules_allT-10000MU/"


labels=["Tx","Ty","Tz"]
positions=["x_global","y_global","z_global"]
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

align_output=open_alignment(thisfile)

plotOn=True

class scifi_layer:
    x_bins, y_bins = 10, 2
    x_start, x_end = -25, 25
    y_start, y_end = -20, 20
    extent = [x_start, x_end, y_start, y_end]

    def __init__(self):
        self.cmax = self.cmin = None
        self.colormap = plt.get_cmap('inferno')
        self.textcolor = 'teal'
        self.p_modules = []
        self.p_mats = []
        # Reminder of quarter arrangement
        #  3  |   2
        # ------------
        #  1  |   0

        for q in [0, 1, 2, 3]:  # quarters
            for m in [0, 1, 2, 3, 4]:  # modules
                self.p_modules.append(patches.Rectangle(
                    (
                        (scifi_layer.x_start if q in [
                         0, 3] else 0) + (scifi_layer.x_bins/2.)*m,
                        (0. if q in [0, 1] else scifi_layer.y_start)
                    ), scifi_layer.x_bins/2., scifi_layer.y_end, linewidth=3, edgecolor='grey', facecolor='none'))
            for mat in range(4*5):
                self.p_mats.append(patches.Rectangle(
                    (
                        (scifi_layer.x_start if q in [
                         0, 3] else 0) + (scifi_layer.x_bins*mat)/4.,
                        (0. if q in [0, 1] else scifi_layer.y_start)
                    ), (scifi_layer.x_bins/2.)/4., scifi_layer.y_end, linewidth=1, edgecolor='grey', facecolor='none', linestyle='dashed'))

        self.map_alignable_to_canvas = {}
        for q, qname in enumerate(["Quarter0", "Quarter1", "Quarter2", "Quarter3"]):  # quarters
            # modules
            for m, mname in enumerate(["Module0", "Module1", "Module2", "Module3", "Module4"]):
                self.map_alignable_to_canvas[(
                    0 if q in [2, 3] else 1,
                    4 - m if q in [1, 3] else m + 5
                )] = qname + mname
                for mat, matname in enumerate(["Mat0", "Mat1", "Mat2", "Mat3"]):  # mats
                    self.map_alignable_to_canvas[(
                        0 if q in [2, 3] else 1,
                        # 20 mats per quadrant
                        (4-m)*4 + mat if q in [1, 3] else 20 + m*4 + mat, -99)] = qname + mname + matname

        self.map_alignable_to_xzposition = {}
        for s, sname in enumerate(["T1", "T2", "T3"]):
            for l, lname in enumerate(["LayerX1", "LayerU", "LayerV", "LayerX2"]):
                for q, qname in enumerate(["Quarter0", "Quarter1", "Quarter2", "Quarter3"]):
                    for m, mname in enumerate(["Module0", "Module1", "Module2", "Module3", "Module4"]):
                        for mat, matname in enumerate(["Mat0", "Mat1", "Mat2", "Mat3"]):
                            # for now vertical modules end up to the same point
                            myx = (3.2 + mat*6.5) + m * \
                                26 if q in [0, 3] else 130 - \
                                m*26 - (3.2 + mat*6.5)
                            myz = 13 + (34.1*s + 3*l)
                            self.map_alignable_to_xzposition[sname +
                                                             lname+qname+mname+matname] = (myz, myx)

    def draw_modules(self, ax):
        for m in self.p_modules:
            new_m = copy.copy(m)
            ax.add_patch(new_m)

    def draw_mats(self, ax):
        for m in self.p_mats:
            new_m = copy.copy(m)
            ax.add_patch(new_m)

    # modded to accept python dictionary
    def plot_data_modules(self, data, title, filename, xlabel='x [a.u.]', ylabel='y [a.u.]', zlabel='',layerinfo=None, save=True):
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.subplots_adjust(left=0.04, bottom=0.04)
#        ax.autoscale()

        matrix_data = np.zeros((2, 10))
        jump_x = (scifi_layer.x_end - scifi_layer.x_start) / \
            (2.*scifi_layer.x_bins)
        jump_y = (scifi_layer.y_end - scifi_layer.y_start) / \
            (2.*scifi_layer.y_bins)
        x_positions = np.linspace(
            start=scifi_layer.x_start, stop=scifi_layer.x_end, num=scifi_layer.x_bins, endpoint=False)
        y_positions = np.linspace(
            start=scifi_layer.y_start, stop=scifi_layer.y_end, num=scifi_layer.y_bins, endpoint=False)
        for y_idx, yv in enumerate(y_positions):
            for x_idx, xv in enumerate(x_positions):
                name_alignables = self.map_alignable_to_canvas[(y_idx, x_idx)]
                if layerinfo:
                    label = data[f'{layerinfo}{name_alignables}']
                else:
                    label=data[f'{name_alignables}']
                matrix_data[y_idx, x_idx] = label
                text_x = xv + jump_x
                text_y = yv + jump_y
                ax.text(text_x, text_y, "%3.2f" % label, ha='center', va='center',
                        fontsize=24, color=self.textcolor, fontweight='bold', rotation=75)
        im = ax.imshow(matrix_data, cmap=self.colormap, origin='top',  # viridis
                       extent=scifi_layer.extent,
                       vmin=self.cmin, vmax=self.cmax,
                       interpolation=None)
        cb = plt.colorbar(im, fraction=0.15, pad=0.02)
        cb.set_label(zlabel)
        ax.set_title(title, fontdict={"fontsize": 24})
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        self.draw_modules(ax)
        ax.set_axis_off()
        
        if save:
            fig.savefig(filename)
            ax.cla()
            fig.clf()
            del fig
            del ax
        else:
            return fig

    def plot_data_mats(self, data, title, filename, xlabel='x [a.u.]', ylabel='y [a.u.]', zlabel=''):
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.subplots_adjust(left=0.04, bottom=0.04)
#        ax.autoscale()

        matrix_data = np.zeros((2, 10*4))
        jump_x = (scifi_layer.x_end - scifi_layer.x_start) / \
            (2.*scifi_layer.x_bins*4)
        jump_y = (scifi_layer.y_end - scifi_layer.y_start) / \
            (2.*scifi_layer.y_bins)
        x_positions = np.linspace(
            start=scifi_layer.x_start, stop=scifi_layer.x_end, num=scifi_layer.x_bins*4, endpoint=False)
        y_positions = np.linspace(
            start=scifi_layer.y_start, stop=scifi_layer.y_end, num=scifi_layer.y_bins, endpoint=False)
        increment = 1
        for y_idx, yv in enumerate(y_positions):
            for x_idx, xv in enumerate(x_positions):
                name_alignables = self.map_alignable_to_canvas[(
                    y_idx, x_idx, -99)]
                label = data[data.index.str.match(r'.+'+name_alignables+'.+')]
                matrix_data[y_idx, x_idx] = label
                text_x = xv + jump_x
                text_y = yv + jump_y/4.*increment + jump_y*0.5
                increment += 1
                if increment > 4:
                    increment = 1
                ax.text(text_x, text_y, "%3.2f" % label, ha='center', va='center',
                        fontsize=18, color=self.textcolor, fontweight='bold', rotation=90)
        im = ax.imshow(matrix_data, cmap=self.colormap, origin='top',  # viridis
                       extent=scifi_layer.extent,
                       vmin=self.cmin, vmax=self.cmax,
                       interpolation=None)
        cb = plt.colorbar(im, fraction=0.15, pad=0.02)
        cb.set_label(zlabel)
        ax.set_title(title, fontdict={"fontsize": 24})
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        self.draw_modules(ax)
        self.draw_mats(ax)
        ax.set_axis_off()
        fig.savefig(filename)
        ax.cla()
        fig.clf()
        plt.close(fig)
        del fig
        del ax

print("###### plot generation ######")
print("..... wait .....")
SF = scifi_layer()

#for non-pandas build: need to remap format of dictionary
#alignables=align_output.keys()
#layertest={}
#layercrop="T1LayerX1"
#labels=["Tx"]
#for label in labels:
#    layertest[label]={}
#    for alignable in alignables:
#        if layercrop in alignable:
#            quarters=alignable.split("|")
#            module=alignable.split("Module")
            
#            print('q', quarters)
#            print('m', module)
#            print('q[0]', quarters[0])
#            print('m[1]', module[1])
            #split quarters for drawing
#            newalignablefirst=f"Quarter{quarters[0][-1]}Module{module[1][0]}"
#            newalignablesecond=f"Quarter{quarters[1][0]}Module{module[1][0]}"
#            layertest[label][newalignablefirst]=1000*align_output[alignable][label][-1]
#            layertest[label][newalignablesecond]=1000*align_output[alignable][label][-1]

#print(layertest["Tx"])

#values=layertest["Tx"].values()
#values=list(values)

#mods_max = max( abs(max(values)), abs(min(values)) )
#SF.textcolor = 'black'
#SF.cmax = mods_max
#SF.cmin = -mods_max
#SF.colormap = plt.get_cmap('RdBu_r')
#result=SF.plot_data_modules(layertest["Tx"], "Modules in T1LayerX1",
#                     f'{outdir}/T1LayerX1Modules-Tx-grid.pdf', zlabel=r'Tx $[\mu m]$',save=False)

#plt.show()

############################
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
        #    3  |   2
        # A ------------ C
        #    1  |   0

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
            im = ax.imshow(matrix_data, cmap=self.colormap, origin='top',  # viridis
                           extent=scifi_all.extent,
                           vmin=self.cmin, vmax=self.cmax,
                           interpolation=None)
            ax.set_title(f'{title},{layer}' , fontdict={"fontsize": 12})
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            self.draw_modules(ax)

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
                            label = data[f'{station}{layer}{name_alignables}']
                        matrix_data[y_idx, x_idx] = label

                        text_x = xv + jump_x
                        text_y = yv + jump_y
                        ax.text(text_x, text_y, "%3.2f" % label if label is not None else "", ha='center', va='center',
                                fontsize=12, color=self.textcolor, fontweight='bold', rotation=75)
                im = ax.imshow(matrix_data, cmap=self.colormap, origin='top',  
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

print("###### plot generation ######")
print("..... wait .....")
SF = scifi_all()

#for non-pandas build: need to remap format of dictionary
#thisfile="20210412_Modules_allT-10000MU/AlignmentResults/parsed_log.json"
align_output=open_alignment(thisfile)
alignables=align_output.keys()

# Separate quarters
def redefine_alignable_format(alignables):
    alignables=[alignable for alignable in alignables if "Module" in alignable]

    layertest={}
    labels=["Tx","Tz"]
    for label in labels:
        layertest[label]={}
        for alignable in alignables:
            intro=alignable.split("Quarter")[0]
            quarters=alignable.split("|")
            module=alignable.split("Module")
            
            #split quarters for drawing
            newalignablefirst=f"{intro}Quarter{quarters[0][-1]}Module{module[1][0]}"
            newalignablesecond=f"{intro}Quarter{quarters[1][0]}Module{module[1][0]}"
            
            layertest[label][newalignablefirst]=1000*align_output[alignable][label][-1]
            layertest[label][newalignablesecond]=1000*align_output[alignable][label][-1]
    return layertest

layertestMU=redefine_alignable_format(alignables)

#thisfile="20210412_Modules_allT-10000MD/AlignmentResults/parsed_log.json"
align_output=open_alignment(thisfile)
alignables=align_output.keys()
layertestMD=redefine_alignable_format(alignables)

#print(layertestMU["Tx"])
print(layertestMU)
values=list(layertestMU["Tx"].values())+list(layertestMD["Tx"].values())
mods_max = max( abs(max(values)), abs(min(values)) )
#SF.cmax = mods_max
#SF.cmin = -mods_max

#manual setting
SF.cmax = 40
SF.cmin=-40

SF.textcolor = 'black'
SF.colormap = plt.get_cmap('BrBG_r')
result=SF.plot_data_modules_multistation(layertestMU["Tx"], "Tx",
                     f'MUModules-Tx-grid.pdf', zlabel=r'Tx $[\mu m]$',save=False)

print("###### plot generation ######")
print("..... wait .....")
SF = scifi_all()

#for non-pandas build: need to remap format of dictionary
#thisfile="20210505_StationsLayers_TxTz_v17r0-3000MU/AlignmentResults/parsed_log.json"
align_output=open_alignment(thisfile)
alignables=align_output.keys()

# Separate quarters
def redefine_alignable_format(alignables):
    alignables=[alignable for alignable in alignables if "Side" in alignable]
    
    layertest={}
    labels=["Tx","Tz"]
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
                
    return layertest

layertest=redefine_alignable_format(alignables)

#print(layertest["Tx"])
values=list(layertest["Tx"].values())
mods_max = max( abs(max(values)), abs(min(values)) )
SF.cmax = mods_max
SF.cmin = -mods_max

#manual setting
#SF.cmax = 40
#SF.cmin=-40

SF.textcolor = 'black'
SF.colormap = plt.get_cmap('BrBG_r')
result=SF.plot_data_modules_multistation(layertest["Tx"], "Tx",
                     f'MUModules-Tx-grid.pdf', zlabel=r'Tx $[\mu m]$',save=False)


values=list(layertest["Tz"].values())
mods_max = max( abs(max(values)), abs(min(values)) )
SF.cmax = mods_max
SF.cmin = -mods_max

result=SF.plot_data_modules_multistation(layertest["Tz"], "Tz",
                     f'MUModules-Tz-grid.pdf', zlabel=r'Tz $[\mu m]$',save=False)


plt.show()
