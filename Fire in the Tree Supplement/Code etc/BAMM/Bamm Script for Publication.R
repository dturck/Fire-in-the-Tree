require(BAMMtools)
require(phytools)

getwd()

tree <- read.tree("lesfinally2018.tre")

plotTree(tree, fsize=0.1, lwd=0.1)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree),
           frame="none",adj=c(1.1,-0.4), cex=0.1)

#event txt files event_data_man6.txt, event_data_man7.txt, event_data_man8.txt

events <-read.csv("event_data_man6.txt")
head(events)

#mcmc output files mcmc_out_man6.txt, mcmc_out_man7.txt, mcmc_out_man8.txt

mcmcout <- read.csv("mcmc_out_man6.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)


#discard 25% of burnin
ed <- getEventData(tree, events, burnin = 0.25)
head(ed$eventData)

shift_probs <- summary(ed)

#?plot.bammdata
bamm.conif <- plot.bammdata(ed, lwd=1, labels =T, cex = 0.1)
addBAMMshifts(ed, cex = 0.3)
addBAMMlegend(bamm.conif)

msc.set <- maximumShiftCredibility(ed, maximize='product')
msc.config <- subsetEventData(edata, index = msc.set$sampleindex)
plot.bammdata(msc.config, lwd=2)
addBAMMshifts(msc.config, cex = 2)

best <- getBestShiftConfiguration(ed, expectedNumberOfShifts=1)
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex=2.5)
#get clade rates, not obvious clades in comments

Pinus_sub_Pinus <- getCladeRates(ed, node=419)
mean(Pinus_sub_Pinus$lambda)

#Ponderosae

Pinus_pon <- getCladeRates(ed, node=448)
mean(Pinus_pon$lambda)


Abies <- getCladeRates(ed, node=349)
mean(Abies$lambda) 

#Hesperocyparis

Hespero <- getCladeRates(ed, node=17)
mean(Hespero$lambda)

Juniperus<- getCladeRates(ed, node=42)
mean(Juniperus$lambda)

Callitriodeae <- getCladeRates(ed, node=108)
mean(Callitriodeae$lambda)

Pinaceae <- getCladeRates(ed, node=344)
mean(Pinaceae$lambda)

#South American Araucariaceae

AraucariaceaeSA <- getCladeRates(ed, node=325)
mean(AraucariaceaeSA$lambda)

#Australasian Araucariaceae

AraucariaceaeAUS <- getCladeRates(ed, node=313)
mean(AraucariaceaeAUS$lambda)

Podocarpaceae <- getCladeRates(ed, node=167)
mean(Podocarpaceae$lambda)

Podocarpus <- getCladeRates(ed, node=174)
mean(Podocarpus$lambda)

#Outgroup to Podocarpus

PodoOut <- getCladeRates(ed, node=287)
mean(PodoOut$lambda)

Taxaceae <- getCladeRates(ed, node=138)
mean(Taxaceae$lambda)

Cupressaceae <- getCladeRates(ed, node=5)
mean(Cupressaceae$lambda)

#Larix and Pseudotsuga

Larpsudotsug <- getCladeRates(ed, node=568)
mean(Larpsudotsug$lambda)

#Chamecyparis

Chamecyp <- getCladeRates(ed, node=96)
mean(Chamecyp$lambda)

#Pinus subg. strobus

Strobus <- getCladeRates(ed, node=491)
mean(Strobus$lambda)

#mcmcout <- read.csv("mcmc_out_man.txt", header=T)
#plot(mcmcout$logLik ~ mcmcout$generation)
#
#burnstart <- floor(0.1 * nrow(mcmcout))
#postburn <- mcmcout[burnstart:nrow(mcmcout), ]
#
#library(coda)
#effectiveSize(postburn$N_shifts)
#effectiveSize(postburn$logLik)
#
#
#edata <- getEventData(tree, events, burnin=0.1)
#cmat <- getCohortMatrix(edata)
#cohorts(cmat, edata)
#
#
#edata <- getEventData(tree, events, burnin=0.1)
#css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=1)
#plot(css)
#
#css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set#.limit = 0.95)
#
#summary(css)
#
#css$number.distinct
#
#plot.credibleshiftset(css)
#
#
#
##rate through time all
#plotRateThroughTime(ed, intervalCol = "red", avgCol = "red", ylim =c(0,1), #cex.axis=2)
#text(x=30, y =0.8, labels = "All conifers", font=4, cex = 2.0, pos=4)
#
#plotRateThroughTime(ed, intervalCol="blue", avgCol="blue", ylim=c(0,1),cex#.axis=1.5)
#text(x=30, y= 0.8, label="Dolphins only", font=4, cex=2.0, pos=4)
#
#?plotRateThroughTime
#
#plotRateThroughTime(ed, intervalCol = "blue", avgCol = "blue", ylim =c(0,1#), cex.axis=2)
#text(x=30, y =0.8, labels = "All Whales", font=4, cex = 2.0, pos=4)
#
#
#tip.rates <- getTipRates(ed)
#str(tip.rates)
#hist(tip.rates$lambda.avg, xlab = "avlambda")
#hist(tip.rates$mu.avg, xlab = "avmu")
#
#
#
