# AI based data fusion model for cancer subtype diagnosis
## Background
The code contained in this repository was assembled and created for the Bachelor Thesis of Floor Laagland. The thesis was on cancer subtype diagnosis using AI based data fusion model for multi-omics data integration. Portion of the code are borrowed from previous work. The implementation of MDICC was used, with slight alterations, according to Yang et al. [1], see `demo.R`, `NetworkFusion.R`, `label.py`, `score.py` and `LocalAffinityMatrix.py`. The implementation of Dominant Set clustering was used, with slight alterations, according to Vascon et al. [2], see `demo_mdicc.m`. The data that is used is available on The Cancer Genome Atlas (TCGA), we performed this analysis on Kidney Renall cell Carcinoma (KIRC) and Liver Hepatocellular Carcinoma (LIHC).

## Practical usage
The steps below show how to run our final model for one cancer type. The final model includes feature selection methods and Dominant Set clustering. For TYPE in file names add the cancer type and for FS add the feature selection, either fsvt for VarianceThreshold or fsk for SelectKBest.

Step 1: Download the data sets for a cancer type.

Step 2: Load these to `feature_selection_TYPE.ipynb` in Jupyter Notebook, using Python version 3.8.12, by changing the paths to your directory.

Step 3: Change the paths to your directory where you want to save the data sets after feature selection.

Step 4: Set up the virtual environment to use MDICC, see https://github.com/yushanqiu/MDICC.

Step 5: Change all paths to your directories in `demo.R`.

Step 6: Make sure that the following line of code in the read data part is commented out, when feature selection is performed, before running the code in RStudio.

```
data <- data[-1]
```

Step 7: Change the directory where the files are saved to your directory.

Step 8: Load `aff_test_S` in MATLAB.

```
A = readmatrix('your_path\afffus_TYPE_FS.csv'); 
```

Step 9: Load `full_clust_TYPE_FS.txt` and `label.csv`to `evaluation_clustering.py` in PyCharm. Specify for TYPE the cancer type and for FS the feature selection method.

```
true_labels = pd.read_csv('your_path/label.csv')

test_label = np.loadtxt("your_path/full_clust_TYPE_FS.txt")[:, :]
```

Step 10: Load `ari_scores_TYPE_FS.csv` and `ari_scores_TYPE_FS.csv` to `comparisons\comparison4.m`.

## References
[1] Y. Yang, S. Tian, Y. Qiu, P. Zhao, and Q. Zou, “Mdicc: Novel method for multi-
omics data integration and cancer subtype identification,” Briefings in Bioinformat-
ics, vol. 23, no. 3, p. bbac132, 2022. https://github.com/yushanqiu/MDICC

[2] S. Vascon, S. R. Buló, V. Murino, and M. Pelillo, “Dslib: An open source library for
the dominant set clustering method,” 2020. https://github.com/xwasco/DominantSetLibrary
