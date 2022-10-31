import os

output_list = open('uniform.cache', 'w')
for f in os.listdir('/global/cscratch1/sd/bkingast/EOS_inference/EOS/LANL_Project_eos/uniform'):
    print(os.path.join('/global/cscratch1/sd/bkingast/EOS_inference/EOS/LANL_Project_eos/uniform', f), end="\n", file=output_list)

output_list.close()

