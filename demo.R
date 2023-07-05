source('C:/Users/floor/OneDrive/Documents/MDICC/MDICC-main/MDICC-main/NetworkFusion.R')
library(Rcpp)
library(parallel)
library(Matrix)

setwd("C:/Users/floor/OneDrive/Documents/MDICC/MDICC-main/MDICC-main")
dyn.load("projsplx_R.dll")

# connect anaconda environment
library(reticulate)
use_virtualenv("C:\\Users\\floor\\OneDrive\\Documents\\MDICC\\MDICC-main\\MDICC-main\\MDICCenv")
Sys.setenv(RETICULATE_PYTHON="C:\\Users\\floor\\OneDrive\\Documents\\MDICC\\MDICC-main\\MDICC-main\\MDICCenv\\Scripts\\python.exe", required=TRUE)
use_python("C:\\Users\\floor\\OneDrive\\Documents\\MDICC\\MDICC-main\\MDICC-main\\MDICCenv\\Scripts\\python.exe", required=TRUE)
py_config()
py_available()
setwd("C:/Users/floor/OneDrive/Documents/MDICC/MDICC-main/MDICC-main")

# connect Python scripts
source_python("LocalAffinityMatrix.py")
source_python("score.py")
source_python("label.py")

# read data
setwd("C:/Users/floor/OneDrive/Documents/MDICC/data/KIRC/data_fs") # data path
list <- list.files()
data <- data.frame()
data1 <- list()
X <- list()
for(i in list){
  path <- i
  data <- read.csv(file = path, header = TRUE)
  data <- data[-1] # comment out this line, if feature selection has been performed
  data11 <- as.matrix(data)
  data1[[i]] = scale(data11, center=TRUE, scale=TRUE) # scale data
  data2 = t(data1[[i]])
  d1 = dist(data2)
  d1 = as.matrix(d1)
  X[[i]] <- d1
}

# parameter setting
k1 = 18 # the neighbor of affinity matrix
k2 = 42 # number of nearest neighbors of hyper parameters of gamma
k3 = 2  # number of cluster
c  = 3  # c = k3(c>2) or c = 3(c=2)

# calculate affinity matrix
aff = list()
for(i in 1:3){
  a = as.matrix(X[[i]])
  xxx = testaff(a,k1)
  aff[[i]] = xxx
}


# network fusion
test = MDICC(aff,c = c,k = k2)
test_S = as.matrix(test)  

# calculate affinity matrix of fused network
aff_test_S = testaff(test_S, k1) # used for DS clustering
write.csv(aff_test_S,file="C:/Users/floor/OneDrive/Documents/MDICC/MDICC-main/MDICC-main/kirc/afffus_kirc_fsk.csv")

# result
score = MDICCscore(test_S,k3,'C:/Users/floor/OneDrive/Documents/MDICC/data/KIRC/label.csv', 'class1') # label path
names(score) = c('RI','ARI','NMI','Accu','F1')
label = MDICClabel(test_S,k3)
MDICCresult = list()
MDICCresult[['score']] = score
MDICCresult[['label']] = label
MDICCresult[['score']] # show scores
MDICCresult[['label']] # show labels
write.csv(MDICCresult[['label']],file="C:/Users/floor/OneDrive/Documents/MDICC/MDICC-main/MDICC-main/kirc/featureselection/labels_km_kirc_k.csv")
write.csv(MDICCresult[['score']],file="C:/Users/floor/OneDrive/Documents/MDICC/MDICC-main/MDICC-main/kirc/featureselection/scores_km_kirc_k.csv")
