import re
import json
import argparse

# filenames must be given in command line
#filename_in = "~/output_alignment/MD_events/20210511/AlignmentResults/alignlog.txt"
#filename_out = "~/output_alignment/MD_events/20210511/Modules_allT_1000MD/parsedlog.json"
#filename_out_extra = "~/output_alignment/MD_events/20210511/Modules_allT_1000MD/globalpars.txt"

parser = argparse.ArgumentParser()
parser.add_argument("filename_in")
parser.add_argument("filename_out")
parser.add_argument("filename_out_extra")
args=parser.parse_args()

filename_in=args.filename_in
filename_out=args.filename_out
filename_out_extra=args.filename_out_extra

print(f"Reading: {filename_in}")

### regex per iteration
regex_iteration=re.compile("Iteration*")
regex_leadingline=re.compile("Values of constraint equation*")
regex_convergence=re.compile("Convergence*")

### regex per alignable
regex_alignable=re.compile("Alignable*")
regex_position_global=re.compile("Global position*")
regex_tracks=re.compile("Number of tracks/hits/outliers*")
regex_alignpars=re.compile("align pars*")
# total chi2 / dof
regex_chi2Dof=re.compile("Total chisquare/dofs*")
# total local ch2 / dofs
#regex_localChi2Dof=re.compile("total local delta chi2 / dof*")

alignments={}
alignments["converged"]=[]
alignments["total_chi2_vals"]=[]  # seperate chi2 value from dfs value!
alignments["total_chi2_dofs"]=[]
labels=["Tx","Ty","Tz","Rx","Ry","Rz"]
with open(filename_in,"r") as inputfile:
    with open(filename_out_extra, "w") as globalparsfile:
        thisIter=-1
        writeGlobal=0
        thisObject=""
        for line in inputfile:
            # properties per iteration
            if regex_iteration.search(line):
                thisIter+=1
                continue
            if regex_leadingline.search(line):
                writeGlobal=1
                globalparsfile.write(f"Iteration: {thisIter}\n")
                globalparsfile.write(line)
                continue
            if regex_convergence.search(line):
                writeGlobal=0
                if "Not" in line:
                    alignments["converged"].append(0)
                    print(f"Notice: alignment not converged in iteration {thisIter}")
                    globalparsfile.write(f"Notice: alignment not converged in iteration {thisIter}")
                    continue
                else:
                    alignments["converged"].append(1)
                    continue
            if regex_chi2Dof.search(line):
                textlist=re.split(r"[\s;/]",line)
                splitted=' '.join(textlist).split()
                name=' '.join(splitted[0:3])
                print("textlist: ", name, splitted[3], "/", splitted[4])
                alignments["total_chi2_vals"].append(float(splitted[3]))
                alignments["total_chi2_dofs"].append(float(splitted[4]))
                continue
            if writeGlobal:
                globalparsfile.write(line)

            # properties per alignable
            if regex_alignable.search(line):
                text,thisObject=line.split(":")
                thisObject=thisObject.strip()
                if thisObject in alignments.keys():
                    continue
                else:
                    alignments[thisObject]={label:[] for label in labels}
                    alignments[thisObject]["x_global"]=[]
                    alignments[thisObject]["y_global"]=[]
                    alignments[thisObject]["z_global"]=[]
                    alignments[thisObject]["nTracks"]=[]
                    alignments[thisObject]["nHits"]=[]
#                    alignments[thisObject]["total_chi2_dofs"]=[]
#                    alignments[thisObject]["local_chi2_dofs"]=[]
                    continue
            if regex_position_global.search(line):
                textlist=re.split("\(|,|\)",line)
                alignments[thisObject]["x_global"].append(textlist[1])
                alignments[thisObject]["y_global"].append(textlist[2])
                alignments[thisObject]["z_global"].append(textlist[3])
                continue
            if regex_tracks.search(line):
                text,trackvars=line.split(":")
                typeslist=trackvars.split()
                alignments[thisObject]["nTracks"].append(typeslist[0])
                alignments[thisObject]["nHits"].append(typeslist[1])
                continue
#            if regex_chi2Dof.search(line):
#                text,totalchivars=line.split(":")
#                typeslist=totalchivars.split()
#                alignments[thisObject]["total_chi2_dofs"].append(typeslist[0])
#                continue
#            if regex_localChi2Dof.search(line):
#                text,localchivars=line.split(":")
#                typeslist=localchivars.split()
#                alignments[thisObject]["local_chi2_dofs"].append(typeslist[0])
#                continue
            if regex_alignpars.search(line):
                text,alignvars=line.split(":")
                varlist=alignvars.split()
                for (label,alignvar) in zip(labels,varlist):
                    alignments[thisObject][label].append(alignvar)
                continue

f=open(filename_out,"w")
f.write( json.dumps(alignments))
f.close()
