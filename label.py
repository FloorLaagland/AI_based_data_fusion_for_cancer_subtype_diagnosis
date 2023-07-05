from sklearn.cluster import KMeans

def MDICClabel(x,k):
    k = int(k) # number of clusters
    kmeans = KMeans(n_clusters=k,random_state = 0).fit(x) # apply kmeans clustering, .fit(X) to train the model
    py_result = kmeans.predict(x) # make use of the learned label to map and predict the labels for the data to be tested
    test_label = py_result.tolist() # map labels to a list
    return test_label