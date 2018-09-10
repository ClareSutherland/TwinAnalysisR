#--------------initialise--------------
#uses openmx and onyx, uncomment and run to install if needbe
#install.packages("OpenMx")
#install.packages("umx")

require(OpenMx)
require(umx)

#--------------loading data--------------
# open the dataset of Australian twin data
data("twinData")

#--------------set variables--------------
# set manifest as bmi
selDVs = c("bmi")

# set mz and dz twin subsets (note, female)
dz = twinData[twinData$zygosity == "DZFF", ]
mz = twinData[twinData$zygosity == "MZFF", ]

#--------------ACE model--------------
#run ACE Model
ACE = umxACE(selDVs = selDVs, dzData = dz, mzData = mz, sep = "")

#Plot the ACE model (square these numbers to get prop. variance explained)
umxSummary(ACE)

#parameters shows you the free parameters          
parameters(ACE)

#--------------AE model--------------
#run AE model by dropping C
AE = umxModify(ACE, update = "c_r1c1", name = "dropC")

#Plot the AE model (square these numbers to get prop. variance explained)
umxSummary(AE, comparison = ACE)

