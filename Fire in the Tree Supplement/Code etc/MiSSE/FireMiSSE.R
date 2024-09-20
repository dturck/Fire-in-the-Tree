# MiSSE runscript R
# Orlando Schwery 28. Jan. 2023

library(hisse)  # updated to 2.1.10
library(phytools)
library(RColorBrewer)

setwd("C:/Users/oschwery/Documents/NMC-Postdoc/FireAdaptationTurck")  # set working dir

phy <- read.tree("version_14/lesfinally2018.tre")  # load tree

sampf <- length(phy$tip.label)/614  # set sampling fraction based on 614 extant species


# Loop through running MiSSE with 1 to 9 rate categories
max_attempted_rates <- 9
all_MiSSEs <- list()
for (i in 1:max_attempted_rates) {
    turnover <- c(1:i)
    eps <- c(1:i)
    all_MiSSEs[[i]] <- MiSSE(phy, f=sampf, turnover=turnover, eps=eps, dt.threads=parallel::detectCores())
}


# check AIC to see how much rates are improved by additional states
MiSSE_AIC <- c()
MiSSE_AICc <- c()

for (i in 1:length(all_MiSSEs)) {
    MiSSE_AIC <- c(MiSSE_AIC, all_MiSSEs[[i]]$AIC)
    MiSSE_AICc <- c(MiSSE_AICc, all_MiSSEs[[i]]$AICc)
}

# plot to see bend in AIC
plot(MiSSE_AIC)
lines(MiSSE_AIC)
lines(MiSSE_AICc, col='blue')

# calculate delta AIC for convenience
delta_AIC <- MiSSE_AIC-(min(MiSSE_AIC))


# calculate marginal ancestral state reconstruction of 3 rate classes and plot [this plots the estimated net diversification rates, no the reconstructed rate classes, those are below]
three.rate.recon <- MarginReconMiSSE(phy=phy, f=sampf, hidden.states=3, pars=all_MiSSEs[[3]]$solution, n.cores=1, AIC=all_MiSSEs[[3]]$AIC)
plot.misse.states(three.rate.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram", fsize=.25, legend="none")

# same for 4 rate cats, as delta AIC <2 (qualitatively same result)
windows()
four.rate.recon <- MarginReconMiSSE(phy=phy, f=sampf, hidden.states=4, pars=all_MiSSEs[[4]]$solution, n.cores=1, AIC=all_MiSSEs[[4]]$AIC)
plot.misse.states(four.rate.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram", fsize=.25, legend="none")


# Custom function to collate rate estimations from MiSSE runs
get_rates <- function(solution) {
    nstates <- solution$hidden.states
    turnovers <- solution$solution[seq(from=1, by = 2, length.out=nstates)]
    epss <- solution$solution[seq(from=2, by = 2, length.out=nstates)]
    trans_rate <- solution$solution["q0"]
    expected_changes <- trans_rate * sum(solution$phy$edge.length)
    lambda <- turnovers/(1+epss)
    mu <- (epss * turnovers)/(1+epss)
    netdiv <- lambda - mu
    divrates <- rbind(turnovers, epss, lambda, mu, netdiv)
    transitions <- cbind(trans_rate, expected_changes)
    outlist <- list(divrates=divrates, transitions=transitions)
    return(outlist)
}
# get for 3 and 4
rates3 <- get_rates(all_MiSSEs[[3]])
rates4 <- get_rates(all_MiSSEs[[4]])


# Custom function to plot rate estimates of two different sets of categories for easy comparison
rateplots <- function(ratesA, ratesB) {
    par(mfrow=c(3,2))
    for (i in 1:nrow(ratesB$divrates)) {
        plot(sort(ratesB$divrates[i, ]), col="black", pch=20, cex=2, main=rownames(ratesB$divrates)[i], ylab='Rate', xlab="Categories")
        points(sort(c(rep(-1, times=abs(ncol(ratesB$divrates)-ncol(ratesA$divrates))), ratesA$divrates[i, ])), col="blue", pch=20, cex=2)
    }
    plot(c(ratesA$transitions[, "expected_changes"], ratesB$transitions[, "expected_changes"]), col=c("blue", "black"), cex=2, pch=20, xlim=c(0.5, 2.5), ylim=c(0, 1+max(ratesA$transitions, ratesB$transitions)), main="Expcected Changes", ylab="nChanges")
}
# plot for 3 and 4
windows()
rateplots(rates3, rates4)


# Custon function to plot estimated rate categories for all tips and reconstructions for nodes
stateplot <- function(rates, nrates, type) {
    plot(rates$phy, type=type, show.tip.label=FALSE, label.offset=1)
    col <- c("darkcyan", "gold", "maroon2", "darkorchid", "seagreen2")  # the colours need to be in the same order as the trait states
    # tiplabels(pch=22, bg=col[as.numeric(as.factor(discretes[,2]))], cex=2, adj=1)  # we get the tip state coloured based on what the state is in the original data (with some subsetting etc)
    tiplabels(pie = rates$tip.mat[, 2:(nrates+1)], piecol = col, cex = 0.25, adj=1)  # we get the node states coloured based on the reconstruction $lik.anc is the likelihood values for each of the different trait states at 
    nodelabels(pie = rates$node.mat[, 2:(nrates+1)], piecol = col, cex = 0.25)  # we get the node states coloured based on the reconstruction $lik.anc is the likelihood values for each of the different trait states at 
    legend("bottomleft", title="Hidden States", legend=colnames(rates$tip.mat)[2:(nrates+1)], fill=col, horiz=FALSE)
}
# plot for 3 and 4
windows()
stateplot(four.rate.recon, 4, "phylogram")
windows()
stateplot(three.rate.recon, 3, "phylogram")



## save collection of plots to PDF
pdf("version_19/Fire_MiSSE_variable_eps.pdf", width=10, height=30)
# AICs
plot(MiSSE_AIC)
lines(MiSSE_AIC)
lines(MiSSE_AICc, col='blue')
# rates and states for 3 and 4 categories
plot.misse.states(three.rate.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram", fsize=.25, legend="none")
stateplot(three.rate.recon, 3, "phylogram")
plot.misse.states(four.rate.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram", fsize=.25, legend="none")
stateplot(four.rate.recon, 4, "phylogram")
# circular trees
stateplot(three.rate.recon, 3, "fan")
stateplot(four.rate.recon, 4, "fan")
# rate comparisons
rateplots(rates3, rates4)
dev.off()


## Save output files
# AIC table
deAICs <- rbind(MiSSE_AIC, MiSSE_AICc, delta_AIC)
colnames(deAICs) <- eps
write.csv(deAICs, file = "version_19/MiSSE_AICs.csv", row.names=TRUE, col.names=TRUE)
# MiSSE inference output
save(all_MiSSEs, file = "version_19/MiSSE_output.RData")
# Marginal reconstructions
save(three.rate.recon, file = "version_19/MiSSE_3-state_recon.RData")
save(four.rate.recon, file = "version_19/MiSSE_4-state_recon.RData")
