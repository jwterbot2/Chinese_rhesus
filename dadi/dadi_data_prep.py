#################################################################################
#####           dadi data preparation                                       #####
#################################################################################
#!/usr/bin/env python3
"""
Script to prepare the inputs (data dictionary) needed for analyses in dadi, and to generate SFS
"""
import dadi
import matplotlib.pyplot as plt
import sys

datafile = '[/path/to/sequencing/data.]vcf' # chromosome data
popfile = '[/path/to/population/file].txt' # text file that assigns individuals with their population id 
dd = dadi.Misc.make_data_dict_vcf(datafile, popfile)

pop_ids = ["pop1"]
ns_full = [158]  
fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns_full, polarized=False)
fs.to_file('[SFS-file-name].fs')

# plot SFS
fig = plt.figure(figsize=(10, 6))
dadi.Plotting.plot_1d_fs(fs)
plt.title('SFS (folded)')
plt.savefig('[folded_sfs_plot_name].png')
plt.close()
