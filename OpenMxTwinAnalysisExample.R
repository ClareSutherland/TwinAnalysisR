#--------------NOTES--------------
#Clare Sutherland 05.06.18
#following OpenMX UserGuide, page 72-78
#https://openmx.ssri.psu.edu/documentation

#if you think that the optimizer is stopping at a local minimum
#you can try replacing mxRun() with mxTryHard() in your script

#--------------initialise--------------
#uses openmx and onyx, uncomment and run to install if needbe
#install.packages("OpenMx")
#install.packages("onyxR")

#initialise openmx and onyx
#require is designed for use inside other functions
#it returns FALSE and gives a warning if the package does not exist (unlike library, which would error)
require(OpenMx)
require(onyxR)

#ACE model (A = additive genetic, C = common enviromental, E = unique enviromental)
#Classic twin study design, comparing MZ and DZ twins

#--------------loading data--------------
#load and view data (we're using an example twin data set which comes with OpenMx)
data(twinData)
View(twinData)

#--------------set variables--------------
#select manifest variables for analysis (here, BMI)
#NB assumes twin 2 is stored after twin 1!!
#latent aceVars is always the same for an ACE model)
selVars <- c('bmi1','bmi2')
aceVars <- c("A1","C1","E1","A2","C2","E2")

#select MZ and DZ twin data for analysis
#NB just female data - select MZMM and DZMM for male data
mzData <- subset(twinData, zygosity=="MZFF", selVars)
dzData <- subset(twinData, zygosity=="DZFF", selVars)

#--------------set starting values--------------
#generate descriptive statistics for calculation of starting values
#starting value mean: take grand mean of mz and dz T1 and T2 observed phenotypic variables
MeanStartValue <- mean(colMeans(rbind(mzData,dzData), na.rm = T), na.rm = T)

#starting value variation: mean variation in observed phenotypic variables for MZ and DZ T1 and T2, divided by 3, sqrt
varMz <- cov(mzData,use="complete")
varDz <- cov(dzData,use="complete")
VarStartValue <- sqrt((mean(c(varMz[1],varMz[4],varDz[1],varDz[4]), na.rm = T))/3)

#--------------ACE model--------------
#variances of latent (ACE) variables i.e. double headed arrows back to self, constrain to be 1
latVariances <- mxPath(from=aceVars, arrows=2, free=FALSE, values=1)

#means
#means of observed variables
#unidirectional  arrows from triangle to the manifest variable, free but given a plausible start value
#because we use the same label (“mean”) for the two means, they are constrained to be equal
obsMeans <- mxPath(from="one", to=selVars, arrows=1, free=TRUE, values = MeanStartValue, labels="mean")
#means of latent variables
#unidirectional arrows from triangle to each latent variable, fixed at zero
latMeans <- mxPath(from="one", to=aceVars, arrows=1, free=FALSE, values=0)

#paths
#path coefficients for twin 1 (NB assumes twin 2 is stored after twin 1 in selVars list)
pathAceT1 <- mxPath(from=c("A1","C1","E1"), to=selVars[1], arrows=1, free=TRUE, values = VarStartValue, label=c("a","c","e"))
#path coefficients for twin 2 (NB assumes twin 2 is stored after twin 1 in selVars list)
pathAceT2 <- mxPath(from=c("A2","C2","E2"), to=selVars[2], arrows=1, free=TRUE, values = VarStartValue, label=c("a","c","e"))
#covariance between A1 & A2 in MZ twins (NB starting value 1!!!!)
covA1A2_MZ <- mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=1)
#covariance between A1 & A2 in DZ twins (NB starting value 0.5!!!)
covA1A2_DZ <- mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=.5)
#covariance between C1 & C2 (equal enviroments assumption)
covC1C2 <- mxPath(from="C1", to="C2", arrows=2, free=FALSE, values=1)

#Create MZ and DZ models
paths <- list(latVariances, latMeans, obsMeans, pathAceT1, pathAceT2, covC1C2)
modelMZ <- mxModel(model="MZ", type="RAM", manifestVars=selVars, latentVars=aceVars, paths, covA1A2_MZ, mxData(observed=mzData, type="raw"))
modelDZ <- mxModel(model="DZ", type="RAM", manifestVars=selVars, latentVars=aceVars, paths, covA1A2_DZ, mxData(observed=dzData, type="raw"))

#create the -2 loglikelihood function
minus2ll <- mxAlgebra(expression=MZ.fitfunction + DZ.fitfunction,
                                name="minus2loglikelihood")

#Run overall ACE model, having specified -2 loglikelihood and combining Dz and Mz models
obj <- mxFitFunctionAlgebra("minus2loglikelihood")
modelACE <- mxModel(model="ACE", modelMZ, modelDZ, minus2ll, obj)

#Run model
fitACE <- mxRun(modelACE)
sumACE <- summary(fitACE)

#Check model free parameters
sumACE

#--------------ACE model estimates--------------
#Look at various aspects of the model seperately
#Add some labels so the table is understandable
Alabel <- "A"
Clabel <- "C"
Elabel <- "E"
#additive genetic variance, a^2 i.e. heritability
A <- mxEval(a*a, fitACE)
#shared environmental variance, c^2
C <- mxEval(c*c, fitACE)
#unique environmental variance, e^2
E <- mxEval(e*e, fitACE)
#total variance
V <- (A+C+E)
#standardized A (proportion of total variance explained by A)
a2 <- round(A/V,2)
#standardized C (proportion of total variance explained by C)
c2 <- round(C/V,2)
# standardized E (proportion of total variance explained by E)
e2 <- round(E/V,2)

#Table of estimates
#note the second line gives you the proportion of total variance explained by each A,C,E, parameter
estACE <- rbind(cbind("ACE model",Alabel, Clabel, Elabel),cbind("standardised",A,C,E),cbind("prop variance explain",a2,c2,e2)) 
estACE
#-2 log likelihood of ACE model
LL_ACE <- mxEval(fitfunction, fitACE)
LL_ACE

#--------------AE model--------------
#re-specify the paths from C1 to observed variable1 and from C2 to observed variable2 to be fixed to zero
#ie now no common enviroment contribution
modelAE <- mxModel(modelACE, name = "AE")
modelAE <- omxSetParameters(modelAE,labels = "c", free = FALSE,values = 0)

#run and summarise AE model
fitAE <- mxRun(modelAE)
sumAE <- summary(fitAE)
sumAE

#-2 log likelihood of AE model
LL_AE <- mxEval(fitfunction, fitAE)
LL_AE

#comparison likelihood of AE model i.e. does the simple AE model fit worse than the more complex ACE model? 
#difference between two -2 log likelihood models approximates chi2 distribution, thus 
#if the difference between ACE and AE models is around significant, then the ACE model is the preferred model
LRT_ACE_AE <- LL_AE - LL_ACE
LRT_ACE_AE

#get out more detailed info from AE model
Alabel <- "A"
Clabel <- "C"
Elabel <- "E"
A_AEmodel <- mxEval(a*a, fitAE)
C_AEmodel <- mxEval(c*c, fitAE)
E_AEmodel <- mxEval(e*e, fitAE)
V_AEmodel <- (A_AEmodel + C_AEmodel + E_AEmodel)
a2_AE <- round(A_AEmodel/V_AEmodel,2)
c2_AE <- round(C_AEmodel/V_AEmodel,2)
e2_AE <- round(E_AEmodel/V_AEmodel,2)

#note that C now equals 0, which makes sense because we constrained it this way
estAE <- rbind(cbind("AE model",Alabel, Clabel, Elabel),cbind("standardised",A_AEmodel, C_AEmodel, E_AEmodel),cbind("prop variance",a2_AE, c2_AE, e2_AE))
estAE

#--------------Optional: draw figures--------------
#to use this code, uncomment the next two lines
#Note that will need a bit of tweaking to get it to look nice (see powerpoint file or Onyx helpguide)
#onyx(fitACE)
#onyx(fitAE)
