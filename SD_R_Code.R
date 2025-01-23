# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
#                                                                                   					#
# 	R code for:                                                                       				#
#                                                                                   					#
#     	 			SEX DIFFERENCES IN HELPING EFFORT REVEAL THE EFFECT OF FUTURE      				#
#                 		  REPRODUCTION ON COOPERATIVE BEHAVIOUR IN BIRDS			            		#
#                                                                                  					#
#            								Philip Downing, Ashleigh Griffin & Charlie Cornwallis	#
#                                                                                   					#
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------- #
									 ### PACKAGES AND FUNCTIONS ###									
# ------------------------------------------------------------------------------------------------- #

# packages

library(ape)
library(coda)
library(MCMCglmm)
library(QuantPsyc)
library(doBy)
library(orthopolynom)
library(metafor)
library(caper)
library(car)

# functions from Nakagawa et al (2015) Meta-analysis of variation: ecological and evolutionary applications and beyond. Methods in Ecol. & Evol. 6, 143-152.
Calc.d <- function(CMean, CVAR, CN, EMean, EVAR, EN, adjusted=T){
                   sPooled <- sqrt(((EN - 1)*EVAR + (CN - 1)*CVAR) / (EN + CN - 2))  
                   d <- (EMean - CMean) / sPooled
                   H.d <- d * (1 - (3 / (4 * (EN + CN - 2) - 1)))
                   if(adjusted==T){return(H.d)}
                   if(adjusted==F){return(d)}}
Calc.SE.d <- function(CN, EN, d){
                      SE <- sqrt(( (CN + EN) / (CN * EN) ) + ( (d^2) / (2 * (CN + EN - 2) ) ) )
                      return(SE)}
# CMean is the mean level of female help
# EMean is the mean level of male help
# CN is the sample size used to calculate the mean level of female help
# EN is the sample size used to calculate the mean level of male help
# CVar is the variance in female help
# EVar is the variance in male help

# other functions
zscore <- function(x){(x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)}
rateDiff <- function(a, b, c, d){rD <- (a/b) - (c/d)
					return(rD)}

# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #

### DATA ###

## Phenotypic Data ##
sexData <- read.csv("...")

## Trees ##
trees <- read.nexus("...")    # contains 20 species
unique(sexData$animal)[which((unique(sexData$animal) %in% trees[[1]]$tip.label) == FALSE)]  # species in data set missing in tree
trees[[1]]$tip.label[which((trees[[1]]$tip.label %in% unique(sexData$animal)) == FALSE)]  # species in tree missing in data set


# ------------------------------------------------------------------------------------------------- #
									  ### DATA MANIPULATION ###										
# ------------------------------------------------------------------------------------------------- #
		
## create effect size ##
sexData$sexDiff <- Calc.d(sexData$femaleX, sexData$var.female, sexData$n.female, sexData$maleX, sexData$var.male, sexData$n.male)
sexData$seD <- Calc.SE.d(sexData$n.female, sexData$n.male, sexData$sexDiff)
sexData$MEV <- sexData$seD ^ 2
sexData$invSE <- 1/sexData$seD

## calculate the rate difference of male and female future breeding in natal group ## 
sexData$inherDiff <- rateDiff(sexData$perMinher, sexData$perMhelp, sexData$perFinher, sexData$perFhelp)

## calculate the rate difference of male and female subordinate breeding in natal group ## 
sexData$Fbred <- 100*(sexData$fOff / sexData$Nchicks)
sexData$Mbred <- 100*(sexData$mOff / sexData$Nchicks)
sexData$reproDiff <- rateDiff(sexData$Mbred, sexData$perMhelp, sexData$Fbred, sexData$perFhelp)

## calculate the percentage of nests with helpers ##
sexData$perCoop <- 100*(sexData$nHelped / (sexData$nHelped + sexData$nNot))

## collapse data set to one entry per species ##
sexDataUnique <- data.frame(summaryBy(sexDiff + MEV + invSE + inherDiff + reproDiff + perCoop + groupSize + perImmigrants + perInher + perFinher + perMinher + Fbred + Mbred ~ animal, keep.names = TRUE, fun = mean, data=sexData))

## calculate weighted mean of mean sex difference in helping effort ##
sexData$N <- sexData$n.female + sexData$n.male
species <- unique(sexData$animal)
numbered <- as.factor(match(sexData$animal, species))
numberedLevels <- levels(numbered)
estimate <- NULL
for(i in 1:length(numberedLevels)) {aSpecies <- which(!is.na(match(numbered, numberedLevels[i])) == TRUE)             
                                    ifelse(length(aSpecies) > 1 , estimate[i] <- (sum(sexData$sexDiff[aSpecies] * sexData$N[aSpecies], na.rm=TRUE)) / (sum(sexData$N[aSpecies], na.rm=TRUE)), estimate[i] <- sexData$sexDiff[aSpecies])}
WM <- data.frame(animal = unique(sexData$animal), estimate = estimate)
sexDataUnique$WM <- WM$estimate[match(sexDataUnique$animal, WM$animal)]

## to get weighted mean of MEV ##
MEVestimate <- NULL
for(i in 1:length(numberedLevels)) {bSpecies <- which(!is.na(match(numbered, numberedLevels[i])) == TRUE)             
ifelse(length(bSpecies) > 1 , MEVestimate[i] <- (sum(sexData$MEV[bSpecies] * sexData$N[bSpecies], na.rm=TRUE)) / (sum(sexData$N[bSpecies], na.rm=TRUE)), MEVestimate[i] <- sexData$MEV[bSpecies])}
WMEV <- data.frame(animal = unique(sexData$animal), MEVestimate = MEVestimate)
sexDataUnique$WMEV <- WMEV$MEVestimate[match(sexDataUnique$animal, WMEV$animal)]


# ------------------------------------------------------------------------------------------------- #
										  ### STATISTICS ###										
# ------------------------------------------------------------------------------------------------- #
																		
#####################################################################################################
################################# sex differences in helping effort #################################
#####################################################################################################

## metafor: mean sex difference in helper investment ##
uniqueModel <- rma(sexDataUnique$WM, vi=sexDataUnique$WMEV)
summary(uniqueModel)    # B = 0.09, lwr = -0.13, upr = 0.31
forest(uniqueModel)

## MCMCglmm: mean sex difference in helper investment accounting for phylogenetic uncertainty ##
priorA <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
MEV <- sexDataUnique$WMEV
INtree <- inverseA(trees[[1]], nodes="TIPS")
interceptModel.start <- MCMCglmm(WM ~ 1, random = ~ animal, data=sexDataUnique, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=priorA, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
interceptModel <- interceptModel.start
for(i in 1:1300){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=interceptModel$Liab[1,], R=interceptModel$VCV[1,3], G=list(G1=interceptModel$VCV[1,1]))
  interceptModel <- MCMCglmm(WM ~ 1, random = ~ animal, data=sexDataUnique, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=priorA, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    interceptModel.start$VCV[i-300,] <- interceptModel$VCV[1,]
    interceptModel.start$Sol[i-300,] <- interceptModel$Sol[1,]
    interceptModel.start$Liab[i-300,] <- interceptModel$Liab[1,]
  }
}
interceptModelA <- interceptModel.start
interceptModelB <- interceptModel.start
interceptModelC <- interceptModel.start
save(interceptModelA, file="...")
save(interceptModelB, file="...")
save(interceptModelC, file="...")
load("...")

# chain convergence
hist(interceptModelA$Liab)
plot(interceptModelA$VCV)     	# animal and units close to 0
plot(interceptModelA$Sol)   	# intercept estimate well mixed
autocorr(interceptModelA$VCV)   # correlation between successive samples < 0.1 for all components
autocorr(interceptModelA$Sol)   # correlation between successive samples < 0.1 for all components
interceptModelSols <- mcmc.list(list(interceptModelA$Sol, interceptModelB$Sol, interceptModelC$Sol))
plot(interceptModelSols)
gelman.diag(interceptModelSols)     # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(interceptModelA$VCV)    # animal failed halfwidth
heidel.diag(interceptModelA$Sol)    # intercept failed halfwidth

# model parameters
summary(interceptModelA)
posterior.mode(interceptModelA$Sol)		# intercept = 0.14
HPDinterval(interceptModelA$Sol)		# intercept = -0.29 to 0.37

# I^2 values (check whether to include MEV)
# phylogenetic heritability
posterior.mode(interceptModelA$VCV[,1] / (interceptModelA$VCV[,1] + interceptModelA$VCV[,2] + interceptModelA$VCV[,3]))
HPDinterval(interceptModelA$VCV[,1] / (interceptModelA$VCV[,1] + interceptModelA$VCV[,2] + interceptModelA$VCV[,3]))
# B = 0.002, lwr = 0.00, upr = 0.19
# units
posterior.mode(interceptModelA$VCV[,3] / (interceptModelA$VCV[,1] + interceptModelA$VCV[,2] + interceptModelA$VCV[,3]))
HPDinterval(interceptModelA$VCV[,3] / (interceptModelA$VCV[,1] + interceptModelA$VCV[,2] + interceptModelA$VCV[,3]))
# B = 0.002, lwr = 0.00, upr = 0.17
# total
posterior.mode((interceptModelA$VCV[,1] + interceptModelA$VCV[,3]) / (interceptModelA$VCV[,1] + interceptModelA$VCV[,2] + interceptModelA$VCV[,3]))
HPDinterval((interceptModelA$VCV[,1] + interceptModelA$VCV[,3]) / (interceptModelA$VCV[,1] + interceptModelA$VCV[,2] + interceptModelA$VCV[,3]))
# B = 0.08, lwr = 0.001, upr = 0.25



#####################################################################################################
########################################## publication bias #########################################
#####################################################################################################

## Trim and Fill ##
TandF <- rma(WM, WMEV, data=sexDataUnique, method="REML")
trimfill(TandF)   # estimated number of missing studies on the left side: 0 (SE = 2.8)
funnel(TandF)


## Egger's test ##
EggerMod <- MCMCglmm(WM ~ invSE, random= ~ animal, data=sexDataUnique, nitt=1100000, burnin=100000, thin=1000, pr=T, prior=priorA, pedigree=trees[[1]], verbose=FALSE)
summary(EggerMod)
posterior.mode(EggerMod$Sol[,1:2])  		# intercept = 0.04; slope = -0.03
HPDinterval(EggerMod$Sol[,1:2])   		# intercept = -0.38 to 0.83; slope = -0.24 to 0.14


#####################################################################################################
########################### helping effort and breeding in the natal group ##########################
#####################################################################################################

## MODEL 1: mean sex difference in helping effort vs. future and subordinate reproduction ##
# remove any NAs
sexDataFull <- sexDataUnique[-which(is.na(sexDataUnique$inherDiff) | is.na(sexDataUnique$reproDiff)),]
# build the model
MEV <- sexDataFull$WMEV
INtree <- inverseA(trees[[1]], nodes="TIPS")
fullModel.start <- MCMCglmm(WM ~ zscore(inherDiff) + zscore(reproDiff), random = ~ animal, data=sexDataFull, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=priorA, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
fullModel <- fullModel.start
for(i in 1:1300){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=fullModel$Liab[1,], R=fullModel$VCV[1,3], G=list(G1=fullModel$VCV[1,1]))
  fullModel <- MCMCglmm(WM ~ zscore(inherDiff) + zscore(reproDiff), random = ~ animal, data=sexDataFull, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=priorA, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    fullModel.start$VCV[i-300,] <- fullModel$VCV[1,]
    fullModel.start$Sol[i-300,] <- fullModel$Sol[1,]
    fullModel.start$Liab[i-300,] <- fullModel$Liab[1,]
  }
}
fullModelA <- fullModel.start
fullModelB <- fullModel.start
fullModelC <- fullModel.start
save(fullModelA, file="...")
save(fullModelB, file="...")
save(fullModelC, file="...")
load("...")

# chain convergence
hist(fullModelA$Liab)
plot(fullModelA$VCV)			# animal and units close to 0
plot(fullModelA$Sol)			# all parameter estimates well mixed
autocorr(fullModelA$VCV)		# correlation between successive samples < 0.1 for all components
autocorr(fullModelA$Sol)		# correlation between successive samples < 0.1 for all components
fullModelSols <- mcmc.list(list(fullModelA$Sol, fullModelB$Sol, fullModelC$Sol))
plot(fullModelSols)
gelman.diag(fullModelSols)		# upper CI = 1.00 for intercept and slope suggesting convergence
heidel.diag(fullModelA$VCV)		# animal and units failed halfwidth
heidel.diag(fullModelA$Sol)		# all parameters passed halfwidth

# model parameters
summary(fullModelA)
posterior.mode(fullModelA$Sol)	# intercept = 0.27, inherDiff = 0.17, reproDiff = 0.30,
HPDinterval(fullModelA$Sol)		# -0.15 to 0.62, -0.25 to 0.46, -0.13 to 0.64


## MODEL 2: mean sex difference in helping effort vs. future reproduction ##
# remove any NAs
sexDataRD <- sexDataUnique[-which(is.na(sexDataUnique$inherDiff)),]
# build the model
MEV <- sexDataRD$WMEV
INtree <- inverseA(trees[[1]], nodes="TIPS")
rdModel.start <- MCMCglmm(WM ~ zscore(inherDiff), random = ~ animal, data=sexDataRD, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=priorA, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
rdModel <- rdModel.start
for(i in 1:1300){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=rdModel$Liab[1,], R=rdModel$VCV[1,3], G=list(G1=rdModel$VCV[1,1]))
  rdModel <- MCMCglmm(WM ~ zscore(inherDiff), random = ~ animal, data=sexDataRD, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=priorA, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    rdModel.start$VCV[i-300,] <- rdModel$VCV[1,]
    rdModel.start$Sol[i-300,] <- rdModel$Sol[1,]
    rdModel.start$Liab[i-300,] <- rdModel$Liab[1,]
  }
}
rdModelA <- rdModel.start
rdModelB <- rdModel.start
rdModelC <- rdModel.start
save(rdModelA, file="...")
save(rdModelB, file="...")
save(rdModelC, file="...")
load("...")

# chain convergence
hist(rdModelA$Liab)
plot(rdModelA$VCV)			# animal and units close to 0
plot(rdModelA$Sol)			# all parameter estimates well mixed
autocorr(rdModelA$VCV)		# correlation between successive samples < 0.1 for all components
autocorr(rdModelA$Sol)		# correlation between successive samples < 0.1 for all components
rdModelSols <- mcmc.list(list(rdModelA$Sol, rdModelB$Sol, rdModelC$Sol))
plot(rdModelSols)
gelman.diag(rdModelSols)			# upper CI = 1.00 for intercept and slope suggesting convergence
heidel.diag(rdModelA$VCV)		# animal failed halfwidth
heidel.diag(rdModelA$Sol)		# intercept and slope passed halfwidth

# model parameters
summary(rdModelA)
posterior.mode(rdModelA$Sol)	# intercept = 0.15, inherDiff = 0.30
HPDinterval(rdModelA$Sol)		# intercept = -0.16 to 0.48, inherDiff = 0.07 to 0.58


## MODEL 3: mean sex difference in helping effort vs. subordinate reproduction ##
sexDataRepro <- sexDataUnique[-which(is.na(sexDataUnique$reproDiff)),]
MEV <- sexDataRepro$WMEV
INtree <- inverseA(trees[[1]], nodes="TIPS")
reproModel.start <- MCMCglmm(WM ~ zscore(reproDiff), random = ~ animal, data=sexDataRepro, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=priorA, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
reproModel <- reproModel.start
for(i in 1:1300){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=reproModel$Liab[1,], R=reproModel$VCV[1,3], G=list(G1=reproModel$VCV[1,1]))
  reproModel <- MCMCglmm(WM ~ zscore(reproDiff), random = ~ animal, data=sexDataRepro, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=priorA, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    reproModel.start$VCV[i-300,] <- reproModel$VCV[1,]
    reproModel.start$Sol[i-300,] <- reproModel$Sol[1,]
    reproModel.start$Liab[i-300,] <- reproModel$Liab[1,]
  }
}
reproModelA <- reproModel.start
reproModelB <- reproModel.start
reproModelC <- reproModel.start
save(reproModelA, file="...")
save(reproModelB, file="...")
save(reproModelC, file="...")
load("...")

# chain convergence
hist(reproModelA$Liab)
plot(reproModelA$VCV)			# ref, animal and units close to 0; species well mixed
plot(reproModelA$Sol)			# all parameter estimates well mixed
autocorr(reproModelA$VCV)		# correlation between successive samples < 0.1 for all components
autocorr(reproModelA$Sol)		# correlation between successive samples < 0.1 for all components
reproModelSols <- mcmc.list(list(reproModelA$Sol, reproModelB$Sol, reproModelC$Sol))
plot(reproModelSols)
gelman.diag(reproModelSols)		# upper CI = 1.00 for intercept and slope suggesting convergence
heidel.diag(reproModelA$VCV)	# animal and units failed halfwidth
heidel.diag(reproModelA$Sol)	# all parameters passed halfwidth

# model parameters
summary(reproModelA)
posterior.mode(reproModelA$Sol)		# intercept = 0.16, reproDiff = 0.34
HPDinterval(reproModelA$Sol)		# intercept = -0.22 to 0.54, reproDiff = 0.07 to 0.62


## MODEL 4: collinearity between fixed effects ##
INtree <- inverseA(trees[[1]], nodes="TIPS")
colModel.start <- MCMCglmm(inherDiff ~ reproDiff, random = ~ animal, data=sexDataFull, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, prior=priorA, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
colModel <- colModel.start
for(i in 1:1300){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
  start <- list(Liab=colModel$Liab[1,], R=colModel$VCV[1,2], G=list(G1=colModel$VCV[1,1]))
  colModel <- MCMCglmm(inherDiff ~ reproDiff, random = ~ animal, data=sexDataFull, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, prior=priorA, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    colModel.start$VCV[i-300,] <- colModel$VCV[1,]
    colModel.start$Sol[i-300,] <- colModel$Sol[1,]
    colModel.start$Liab[i-300,] <- colModel$Liab[1,]
  }
}
colModelA <- colModel.start
colModelB <- colModel.start
colModelC <- colModel.start
save(colModelA, file="...")
save(colModelB, file="...")
save(colModelC, file="...")
load("...")

# chain convergence
hist(colModelA$Liab)
plot(colModelA$VCV)			# animal and units close to 0
plot(colModelA$Sol)			# all parameter estimates well mixed
autocorr(colModelA$VCV)		# correlation between successive samples < 0.1 for all components
autocorr(colModelA$Sol)		# correlation between successive samples < 0.1 for all components
colModelSols <- mcmc.list(list(colModelA$Sol, colModelB$Sol, colModelC$Sol))
plot(colModelSols)
gelman.diag(colModelSols)		# upper CI = 1.00 for intercept and slope suggesting convergence
heidel.diag(colModelA$VCV)		# animal failed halfwidth
heidel.diag(colModelA$Sol)		# all parameters passed halfwidth

# model parameters
summary(colModelA)
posterior.mode(colModelA$Sol)	# intercept = 0.27, slope = 1.55
HPDinterval(colModelA$Sol)		# intercept = -0.01 to 0.63, slope = 0.07 to 3.38


## MODEL 5: female and male routes to breeding ##

sexDataUnique$sdiff <- abs(sexDataUnique$Mbred - sexDataUnique$Fbred)
sexDataUnique$fdiff <- abs(sexDataUnique$perMinher - sexDataUnique$perFinher)

priorB <- list(R = list(V=diag(2), nu=1.002),G=list(G1=list(V=diag(2), nu=1.002)))
INtree <- inverseA(trees[[1]], nodes="TIPS")
fsdModel.start <- MCMCglmm(cbind(sdiff, fdiff) ~ trait - 1, random = ~idh(trait):animal, rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data=sexDataUnique, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, prior=priorB, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
fsdModel <- fsdModel.start
for(i in 1:1300){
  INtree <- inverseA(trees[[i]], nodes="TIPS")
start <- list(Liab=fsdModel$Liab[1,], R=diag(fsdModel$VCV[1,c(1,2)]), G=list(G1=matrix(fsdModel$VCV[1,c(3,4,5,6)], nrow=2)))
  fsdModel <- MCMCglmm(cbind(sdiff, fdiff) ~ trait - 1, random = ~idh(trait):animal, rcov = ~us(trait):units, family = c("gaussian", "gaussian"), data=sexDataUnique, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, prior=priorB, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    fsdModel.start$VCV[i-300,] <- fsdModel$VCV[1,]
    fsdModel.start$Sol[i-300,] <- fsdModel$Sol[1,]
    fsdModel.start$Liab[i-300,] <- fsdModel$Liab[1,]
  }
}
fsdModelA <- fsdModel.start
fsdModelB <- fsdModel.start
fsdModelC <- fsdModel.start
save(fsdModelA, file="...")
save(fsdModelB, file="...")
save(fsdModelC, file="...")
load("...")

# chain convergence
hist(fsdModelA$Liab)
plot(fsdModelA$VCV)			# animal and units close to 0
plot(fsdModelA$Sol)			# all parameter estimates well mixed
autocorr(fsdModelA$VCV)		# correlation between successive samples < 0.1 for all components
autocorr(fsdModelA$Sol)		# correlation between successive samples < 0.1 for all components
fsdModelSols <- mcmc.list(list(fsdModelA$Sol, fsdModelB$Sol, fsdModelC$Sol))
plot(fsdModelSols)
gelman.diag(fsdModelSols)		# upper CI = 1.00 for intercept and slope suggesting convergence
heidel.diag(fsdModelA$VCV)		# all parameters passed halfwidth
heidel.diag(fsdModelA$Sol)		# all parameters passed halfwidth

# model parameters
summary(fsdModelA)
posterior.mode(fsdModelA$Sol)	# sub repro diff = 4.57, future repro diff = 23.7
HPDinterval(fsdModelA$Sol)		# sub repro diff = -0.64 to 9.12, future repro diff = -2.05 to 53.98


# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #
#                                                                                  					#
#                            			END - thanks for reading! p.      	 	                		#
#                                             					                                    #
# ------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------- #