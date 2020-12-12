# batch renaming of frap files script

import os

thedir = r'Z:\ilya\code\fastGCaMP_analysis\jg8-frap\data\exp2_20201210'

files = os.listdir(thedir)

files_filt = [f for f in files if 'csv' in f]
for f in files_filt:
    f_new = f.replace('.00', '.20201210.00')
    f_new_new= f_new.replace('.stim', '.regular.stim')
    print(f_new_new)
    os.rename(os.path.join(thedir, f), os.path.join(thedir, f_new_new))