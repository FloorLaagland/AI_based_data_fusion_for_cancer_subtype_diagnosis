"""
@author: Floor Laagland [s1054408]
"""

# imports
import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import normalized_mutual_info_score

"""
This Python script computes the ARI and NMI between the true labels and the predicted labels of a cancer type per 
clustering.
"""

if __name__ == '__main__':

    # load true labels
    true_labels = pd.read_csv('C:/Users/floor/OneDrive/Documents/MDICC/data/LIHC/label.csv')
    true_labels = true_labels['class1'].tolist()
    true_labels = np.array(true_labels)

    # load test labels
    test_label = []
    test_label = np.loadtxt("C:/Users/floor/OneDrive/Documents/MDICC/MDICC-main/MDICC-main/lihc/full_clust_lihc_fsvt.txt")[:, :]

    ari_scores = []
    nmi_scores = []
    for i in range(0, len(test_label[1, :])):
        test_label1 = test_label[:, i]

        test_labels = []
        for string in test_label1:
            converted_int = int(string)
            test_labels.append(converted_int)

        # convert values in test_labels to the right value as is used in true_labels
        # for LIHC:
        for i in range(0, len(test_labels)):
            if test_labels[i] == 1:
                test_labels[i] = 1
            elif test_labels[i] >= 2:
                test_labels[i] = 11

        # for KIRC:
        # for i in range(0, len(test_labels)):
        #     if test_labels[i] == 1:
        #         test_labels[i] = 0
        #     elif test_labels[i] >= 2:
        #         test_labels[i] = 1

        pred_labels = test_labels

        # calculate ARI and NMI for the clusterings
        ARI = adjusted_rand_score(pred_labels, true_labels)
        NMI = normalized_mutual_info_score(pred_labels, true_labels)

        # save ARI and NMI in a list
        ari_scores.append(ARI)
        nmi_scores.append(NMI)

    # write scores to .csv file
    np.savetxt("C:/Users/floor/OneDrive/Documents/MDICC/MDICC-main/MDICC-main/kirc/ari_scores_kirc_fsvt.csv", ari_scores, delimiter=", ")
    np.savetxt("C:/Users/floor/OneDrive/Documents/MDICC/MDICC-main/MDICC-main/kirc/nmi_scores_kirc_fsvt.csv", nmi_scores, delimiter=", ")
