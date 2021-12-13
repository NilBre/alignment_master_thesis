### filenames ###
MU_filenames_30k_modules_allDof = []  # usual constraints, TxTzRxRz and also align modules
MU_filenames_config5_10k = []  # config 5 10k events, 2019 data ?
# config 5, 2019 data MU, 5000 events (retest)
MU_filenames_config5_2019 = []
# config 5, 2020 data MU, 10000 events
MU_filenames_config5_2020_fulliter = []
# compare highmomentumtracks and goodlongtracks -> high momentum 
MU_filenames_HE_2019_5000 = []

### dir ###
MU_dir_30k_modules_allDof = "../output/alignment_runs/done_runs/HE_runs/12/iter20/AlignmentResults/"
MU_dir_config5_10k = "../output/data/config5_10k/AlignmentResults/"
MU_dir_config5_2019 = "../output/data/config5_2019/AlignmentResults/"
MU_dir_config5_2020_fulliter = "../output/data/config5_2020_fulliter/AlignmentResults/"  # noch abzuwarten wwie viele iteration es werden
MU_dir_HE_2019_5000 = "../../../output_alignment/20211206/highmomentumttracks/AlignmentResults/"

### concatenate ###
for count in range(10):
  MU_filenames_config5_10k.append(MU_dir_config5_10k + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_config5_2019.append(MU_dir_config5_2019 + f'Iter{count}/alignlog_ft_stationslayers.txt')
  MU_filenames_HE_2019_5000.append(MU_dir_HE_2019_5000 + f'Iter{count}/alignlog_ft_stationslayers.txt')
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
