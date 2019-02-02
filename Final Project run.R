graphics.off()
setwd("D:/study/RMIT/1-2/Applied Bayesian Statistics/Final project")
myData = read.csv("fifa.csv")
yName = "eur_value" ; xName = c("age","eur_wage","overall","international_reputation")
fileNameRoot = "fifa"
numSavedSteps = 10000 ; thinSteps=3

graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Final Project.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
xPred = c(32 ,	565000 ,	94 ,	 5)
#colnames(xPred) = c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]")
mcmcCoda = genMCMC( data=myData , xName=xName , yName=yName , 
                    numSavedSteps=numSavedSteps , thinSteps=thinSteps , 
                    saveName=fileNameRoot , xPred = xPred )
#stopTime = proc.time()
#duration = stopTime - startTime
#show(duration)
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType="png" )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:

summaryInfo = smryMCMC( mcmcCoda , 
                        saveName=fileNameRoot  )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , xName=xName , yName=yName , 
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType="png" )
#------------------------------------------------------------------------------- 
