### filenames ###
MU_filenames_30k_modules_allDof = []  # usual constraints, TxTzRxRz and also align modules
MU_filenames_config5_10k = []  # config 5 10k events, 2019 data ?
# config 5, 2019 data MU, 5000 events (retest)
MU_filenames_config5_2019 = []
# config 5, 2020 data MU, 10000 events
MU_filenames_config5_2020_fulliter = []
# compare highmomentumtracks and goodlongtracks -> high momentum 
MU_filenames_HE_2019_5000 = []

###### alignment runs from stack #######
MU_filenames_config5_2019 = []
MU_filenames_config5_2020_TxRz = []
# stack: test GL tracks for #9 to #12
MU_filenames_no9_10k_2020 = []
MU_filenames_no10_10k_2020 = []
MU_filenames_no11_10k_2020 = []
MU_filenames_no12_10k_2020 = []

# stack: low and normal lumi tests
MU_filenames_low_lumi = []
MU_filenames_normal_lumi = []
MU_filenames_modules_low_lumi = []
MU_filenames_modules_normal_lumi = []

# with clusterbias hack added --> CB for clusterbias
MU_filenames_CB_low_lumi = []
MU_filenames_CB_normal_lumi = []
MU_filenames_CB_modules_low = []
MU_filenames_CB_modules_normal = []

# remove constraints, low lumi
MU_filenames_CB_1st = []
MU_filenames_CB_2nd = []

# misalign + clusterbias + low lumi
a_file = []
b_file = []
c_file = []
d_file = []

### dir ###
MU_dir_30k_modules_allDof = "../output/alignment_runs/done_runs/HE_runs/12/iter20/AlignmentResults/"
MU_dir_config5_10k = "../output/data/config5_10k/AlignmentResults/"
MU_dir_config5_2019 = "../output/data/config5_2019/AlignmentResults/"
MU_dir_config5_2020_fulliter = "../output/data/config5_2020_fulliter/AlignmentResults/"  # noch abzuwarten wwie viele iteration es werden
MU_dir_HE_2019_5000 = "../../../output_alignment/20211206/highmomentumttracks/config5_2019/AlignmentResults/"

### dir from stack ###
MU_dir_config5_2019 = "../output/stack_alignment/config5_2019/AlignmentResults/"
MU_dir_config5_2020_TxRz = "../output/stack_alignment/TxRz_2020_GL/AlignmentResults/"
#
MU_dir_no9_10k_2020 = "../../../../../../../interactive_storage/nbreer/build_stack/stack/Alignment/output/no9to12/9/AlignmentResults/"
MU_dir_no10_10k_2020 = "../../../../../../../interactive_storage/nbreer/build_stack/stack/Alignment/output/no9to12/10/AlignmentResults/"
MU_dir_no11_10k_2020 = "../../../../../../../interactive_storage/nbreer/build_stack/stack/Alignment/output/no9to12/11/AlignmentResults/"
MU_dir_no12_10k_2020 = "../../../../../../../interactive_storage/nbreer/build_stack/stack/Alignment/output/no9to12/12/AlignmentResults/"
#
MU_dir_low_lumi = "../output/stack_alignment/test_low_lumi/AlignmentResults/"
MU_dir_normal_lumi = "../output/stack_alignment/test_normal_lumi/AlignmentResults/"
MU_dir_modules_low_lumi = "../output/stack_alignment/modules_low_lumi/AlignmentResults/"
MU_dir_modules_normal_lumi = "../output/stack_alignment/modules_normal_lumi/AlignmentResults/"

# CB
MU_dir_CB_low_lumi = "../output/stack_alignment/clusterbias_hack/low_lumi/AlignmentResults/"
MU_dir_CB_normal_lumi = "../output/stack_alignment/clusterbias_hack/normal_lumi/AlignmentResults/"
MU_dir_CB_modules_low = "../output/stack_alignment/clusterbias_hack/modules_low/AlignmentResults/"
MU_dir_CB_modules_normal = "../output/stack_alignment/clusterbias_hack/modules_normal/AlignmentResults/"

#
MU_dir_CB_1st = "../output/stack_alignment/clusterbias_hack/without_1st_constr/AlignmentResults/"
MU_dir_CB_2nd = "../output/stack_alignment/clusterbias_hack/without_2nd_constr/AlignmentResults/"

# cb + low + mis
a_dir = "../output/misalignment_runs/clusterbias/0/AlignmentResults/"
b_dir = "../output/misalignment_runs/clusterbias/1/AlignmentResults/"
c_dir = "../output/misalignment_runs/clusterbias/2/AlignmentResults/"
d_dir = "../output/misalignment_runs/clusterbias/3/AlignmentResults/"

### concatenate ###
for count in range(10):
  stl = f'Iter{count}/alignlog_ft_stationslayers.txt'
  MU_filenames_config5_10k.append(MU_dir_config5_10k + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_config5_2019.append(MU_dir_config5_2019 + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_HE_2019_5000.append(MU_dir_HE_2019_5000 + f'Iter{count}/alignlog_ft_stationslayers.txt')
  # for alignment stack
  MU_filenames_config5_2019.append(MU_dir_config5_2019 + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_config5_2020_TxRz.append(MU_dir_config5_2020_TxRz + f'Iter{count}/alignlog_ft_stationslayers.txt')
  #
  MU_filenames_no9_10k_2020.append(MU_dir_no9_10k_2020 + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_no10_10k_2020.append(MU_dir_no10_10k_2020 + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_no11_10k_2020.append(MU_dir_no11_10k_2020 + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_no12_10k_2020.append(MU_dir_no12_10k_2020 + f'Iter{count}/alignlog_ft_stationslayers.txt')
  #
  MU_filenames_low_lumi.append(MU_dir_low_lumi + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_normal_lumi.append(MU_dir_normal_lumi + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_modules_low_lumi.append(MU_dir_modules_low_lumi + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_modules_normal_lumi.append(MU_dir_modules_normal_lumi + f'Iter{count}/alignlog_ft_stationslayers.txt')  
  #
  MU_filenames_CB_low_lumi.append(MU_dir_CB_low_lumi + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_CB_normal_lumi.append(MU_dir_CB_normal_lumi + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_CB_modules_low.append(MU_dir_CB_modules_low + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_CB_modules_normal.append(MU_dir_CB_modules_normal + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_CB_1st.append(MU_dir_CB_1st + stl)
  MU_filenames_CB_2nd.append(MU_dir_CB_2nd + stl)

  a_file.append(a_dir + stl)
  b_file.append(b_dir + stl)
  c_file.append(c_dir + stl)
  d_file.append(d_dir + stl)
for count in range(20):
  MU_filenames_30k_modules_allDof.append(MU_dir_30k_modules_allDof + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_config5_2020_fulliter.append(MU_dir_config5_2020_fulliter + f'Iter{count}/alignlog_ft_stationslayers.txt')

### fill alignlog.txt ###
with open(MU_dir_30k_modules_allDof + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_30k_modules_allDof:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_config5_10k + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_config5_10k:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_config5_2019 + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_config5_2019:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_config5_2020_fulliter + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_config5_2020_fulliter:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_HE_2019_5000 + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_HE_2019_5000:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

# from stack alignment
with open(MU_dir_config5_2019 + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_config5_2019:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_config5_2020_TxRz + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_config5_2020_TxRz:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_no9_10k_2020 + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_no9_10k_2020:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_no10_10k_2020 + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_no10_10k_2020:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_no11_10k_2020 + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_no11_10k_2020:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_no12_10k_2020 + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_no12_10k_2020:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
#######

with open(MU_dir_low_lumi + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_low_lumi:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_normal_lumi + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_normal_lumi:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_modules_low_lumi + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_modules_low_lumi:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_modules_normal_lumi + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_modules_normal_lumi:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

## for clusterbias hack
with open(MU_dir_CB_low_lumi + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_CB_low_lumi:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_CB_normal_lumi + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_CB_normal_lumi:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_CB_modules_low + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_CB_modules_low:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_CB_modules_normal + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_CB_modules_normal:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
###
with open(MU_dir_CB_1st + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_CB_1st:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(MU_dir_CB_2nd + 'alignlog.txt', 'w') as outfile:
    for fname in MU_filenames_CB_2nd:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
########## misalign

with open(a_dir + 'alignlog.txt', 'w') as outfile:
    for fname in a_file:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(b_dir + 'alignlog.txt', 'w') as outfile:
    for fname in b_file:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(c_dir + 'alignlog.txt', 'w') as outfile:
    for fname in c_file:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

with open(d_dir + 'alignlog.txt', 'w') as outfile:
    for fname in d_file:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
