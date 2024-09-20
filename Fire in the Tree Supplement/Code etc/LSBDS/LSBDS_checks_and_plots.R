# Visualise the LSBDS results
# Orlando Schwery, 30. Jan. 2023

library(RevGadgets)
library(coda)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)

setwd("C:/Users/oschwery/Documents/NMC-Postdoc/FireAdaptationTurck")

current_model = "LSBDSext"
current_ncats = "3"


# Read in log files
tree = readTrees("version_14/lesfinally2018.tre")
rates = readTrace(paste("version_14/output_uyedalab2/", current_model, "_", current_ncats, "_rates.log", sep=""))
trace = readTrace(paste("version_14/output_uyedalab2/", current_model, "_", current_ncats, ".log", sep=""))


####### Convergence Assessment (in addition to Tracer)

trace_quant <- readTrace(path = paste("version_14/output_uyedalab2/", current_model, "_", current_ncats, ".log", sep=""), burnin = 0)
trace_quant_MCMC <- as.mcmc(trace_quant[[1]])
ESStrace <- effectiveSize(trace_quant_MCMC)
ESStrace[order(ESStrace)]
traceplot(trace_quant_MCMC)

trace_quant_rates <- readTrace(path = paste("version_14/output_uyedalab2/", current_model, "_", current_ncats, "_rates.log", sep=""), burnin = 0)
trace_quant_MCMC_rates <- as.mcmc(trace_quant_rates[[1]])
ESSrate <- effectiveSize(trace_quant_MCMC_rates)
ESSrate[order(ESSrate)]
traceplot(trace_quant_MCMC_rates)


####### Rate and shift Estimates on Tree

# summarise rate estimates and shifts and plot
combined <- processBranchData(tree    = tree, 
                              dat     = rates,
                              burnin  = 0,
                              summary = "mean",
                              net_div = TRUE)

# make plots for net diversification, speciation, extinction, and rate shifts
netdivtree <- plotTree(combined, color_branch_by = "net_div", branch_color = c("#00008B", "#ffd700"), tip_labels_size = 0.75, lineend="square", line_width=0.9, ladderize=FALSE) #, tree_layout = "circular")
lambdatree <- plotTree(combined, color_branch_by = "avg_lambda", branch_color = c("#00008B", "#ffd700"), tip_labels_size = 0.75, lineend="square", line_width=0.9, ladderize=FALSE) #, tree_layout = "circular")
mutree <- plotTree(combined, color_branch_by = "avg_mu", branch_color = c("#00008B", "#ffd700"), tip_labels_size = 0.75, lineend="square", line_width=0.9, ladderize=FALSE) #, tree_layout = "circular")
shifttree <- plotTree(combined, color_branch_by = "num_shifts", branch_color = c("#00008B", "#ffd700"), tip_labels_size = 0.75, lineend="square", line_width=0.9, ladderize=FALSE) #, tree_layout = "circular")

grid.arrange(netdivtree, shifttree, ncol=2)

#get shift nodes
combined_table <- as.data.frame(combined[[1]][[1]]@data)
combined_table[order(combined_table$num_shifts, decreasing=TRUE), ][1:10, ]

# exploratory plot of posterior shift probabilities
plot(combined_table$num_shifts)
abline(h=c(0.25, 0.5, 0.75), col="blue")


####### Ancestral State Reconstructions of Rate Categories

# read in and process the log file
pdata <- processSSE(paste("version_14/output_uyedalab2/", current_model, "_", current_ncats,".log", sep=""))

# plot the posterior probabilities of rate estimates
ratedensity <- plotMuSSE(pdata)

# read in and process the ancestral states of rate categories
anc_states_file <- paste("version_14/output_uyedalab2/", current_model, "_", current_ncats, "_anc_states.tree", sep="")
p_anc <- processAncStates(path = anc_states_file)

# plot the MAP ancestral states (curcular and regular tree)
statestree_circ <- plotAncStatesMAP(p_anc, tree_layout = "circular", tip_labels_size = 0.75, lineend="square", ladderize=FALSE)

statestree <- plotAncStatesMAP(p_anc, tip_labels_size = 0.75, lineend="square", ladderize=FALSE)

grid.arrange(statestree, shifttree, ncol=2)

# Plot reconstruction of rate categories with probabilities on nodes and tips as pie charts
pietree <- plotAncStatesPie(p_anc,
                        # Add text labels to the tip pie symbols
                        tip_labels_states = TRUE,
                        # Offset those text labels slightly
                        tip_labels_states_offset = .5,
                        # Offset the tip labels to make room for tip pies
                        tip_labels_offset = 1, tip_labels_size = 0.75,
                        # Move tip pies right slightly 
                        tip_pie_nudge_x = .07,
                        # Change the size of node and tip pies  
                        tip_pie_size = 0.8,
                        node_pie_size = 1, 
                        lineend="square", ladderize=FALSE)

grid.arrange(statestree, pietree, shifttree, ncol=3)


#################################
# Save plots to PDF

pdf(file=paste("version_14/results/", current_model, "_", current_ncats, "_ratesTrees.pdf",sep=""), width=20, height=19)
grid.arrange(netdivtree, lambdatree, mutree, shifttree, statestree, ncol=5)
dev.off()

pdf(file=paste("version_14/results/", current_model, "_", current_ncats, "_ratesDens.pdf",sep=""), width=10, height=7)
ratedensity
dev.off()

pdf(file=paste("version_14/results/", current_model, "_", current_ncats, "_catCirc.pdf",sep=""), width=20, height=20)
statestree_circ
dev.off()
