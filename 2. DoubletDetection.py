# Refer to the official course https://doubletdetection.readthedocs.io/en/stable/

import scanpy as sc   # 1.9.5
import numpy as np   # 1.26.0
import doubletdetection   # 4.2
import pandas as pd   # 2.1.1
import os

# change  the current path
os.chdir('E:/Data/GSE220946')
dir_list = os.listdir()

for folder in dir_list:
    print(f"Processing {folder + '.mtx.gz'}...")
    file_path = os.path.join(folder, folder + '.mtx.gz')
    adata = sc.read(file_path)
    raw_counts = np.transpose(adata.X)

    clf = doubletdetection.BoostClassifier()
    labels = clf.fit(raw_counts).predict()
    scores = clf.doublet_score()

    output_file_path = os.path.join(folder, folder + '_doubletdetection.csv')
    arr = scores.mask
    df = pd.DataFrame(arr)
    df.to_csv(output_file_path, index=False)

