import os
import subprocess
# sub_dirs = [x[0] for x in os.walk('/home/CBBI/yangj8/TARGET') if 'log' not in x[0] and 'SRX' not in x[0]][1:]
base_dir = os.path.dirname(os.getcwd())
# bkup_dir = os.path.join(base_dir, 'bkup')
srx_sub = [x[0] for x in os.walk(base_dir) if 'log' not in x[0]][1:]
bb = list()
for i in range(0, len(srx_sub)):
    for file in os.listdir(srx_sub[i]):
        if file.endswith('tar'):
            bb.append(os.path.join(srx_sub[i], file))
            break
for i in range(0, len(bb)):
    for file in os.listdir(bb[i]):
        if file.endswith('tar'):
            print 'moving ' + file + ' to ' + bkup_dir
            command = 'mv ' + os.path.join(bb[i], file) + ' ' + bkup_dir
            p = subprocess.Popen(command, shell=True)
            p.communicate()
for f in bb[1:2]:
    command = 'rm ' + f
    p = subprocess.Popen(command, shell=True)
    p.communicate()
# srx_dir = [x[0] for x in os.walk('/home/CBBI/yangj8/TARGET') if 'SRX' in x[0]][1:]
# for i in range(0, len(sub_dirs)):
#     print 'processing ' + str(i) +'/'+ str(len(sub_dirs))
#     for file in os.listdir(sub_dirs[i]):
#         if file.endswith('.tar'):
#             command = 'tar -xvf ' + os.path.join(sub_dirs[i], file) + ' -C ' + sub_dirs[i]
#             p = subprocess.Popen(command, shell=True)
#             p.communicate()
#     for file in os.listdir(sub_dirs[i]):
#         if file.endswith('.gz'):
#             command = 'gunzip ' + os.path.join(sub_dirs[i], file) + ' < ' + sub_dirs[i]
#             p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
#             p.communicate()
#
# bb = list()
# for i in range(0, len(srx_sub)):
#     for file in os.listdir(srx_sub[i]):
#         if file.endswith('gz'):
#             bb.append(srx_sub[i])
#             break
# aa = list()
# for i in range(0, len(srx_sub)):
#     for file in os.listdir(srx_sub[i]):
#         if file.endswith('fastq'):
#             aa.append(srx_sub[i])
#             break
#
# set(aa).intersection(set(bb))
#
# import os
# only_subs = sub_dirs = [x[0].split('/')[-1] for x in os.walk('/home/CBBI/yangj8/TARGET') if 'log' not in x[0] and 'SRX' not in x[0]][1:]
# import pandas as pd
# all_data = pd.read_csv('new_mani.txt', sep='\t')
# all_data = all_data[~(all_data.id.isin(only_subs))]
# all_data.to_csv('new_mani.txt', sep='\t')
#
#
# bb = list()
# for i in range(0, len(sub_dirs)):
#     print 'processing' + sub_dirs[i]
#     pp = [x[0] for x in os.walk(sub_dirs[i]) if 'SRX' in x[0]]
#     files = [x for x in os.listdir(sub_dirs[i]) if x.endswith('fastq')]
#     if len(pp) > 0 or len(files) > 0:
#         continue
#     else:
#         for file in os.listdir(sub_dirs[i]):
#             if file.endswith('.tar'):
#                 command = 'tar -xvf ' + os.path.join(sub_dirs[i], file) + ' -C ' + sub_dirs[i]
#                 p = subprocess.Popen(command, shell=True)
#                 p.communicate()
#         for file in os.listdir(sub_dirs[i]):
#             if file.endswith('.gz'):
#                 command = 'gunzip ' + os.path.join(sub_dirs[i], file) + ' < ' + sub_dirs[i]
#                 p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
#                 p.communicate()
