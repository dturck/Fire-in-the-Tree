Dataset for: Fire in the Tree: The Origin and Distribution of Fire-adapted Traits within Conifers and their Influence on Speciation Rates across the Conifer Phylogeny.  

Please note .tre files are newick formatted trees from Leslie et al. (2018) An overview of extant conifer evolution from the perspective of the fossil record and Jin et al. (2021) Phylogenomic and ecological analyses reveal the spatiotemporal evolution of global pines.

What everything is and how to use it:

CODE

Ancestral State Reconstruction: 

Leslie et al Ancestral State Recon Files

No root state Conifer Phytools Code for Publication.R – R file to get reconstructions
conif_species_fire_trait_post_res_redo_20_aug_2024.csv – trait file for R script
lesfinally2018.tre – Newick formatted tree from Leslie et al. 2018

Jin et al Ancestral State Recon Files

Jin et al ACE file.R - R file to get reconstructions (basically a clone of the No root state Conifer Phytools Code for Publication.R)
Pinus_fire_for_ace_with_outgroups.csv – trait file for R script
cds_time_16_fossils_newick.tre - Newick formatted tree using 16 fossil calibrations from Jin et al. 2021
cds_time_4_fossils_newick.tre - Newick formatted tree using 16 fossil calibrations from Jin et al. 2021

BAMM

BAMM_Script_divcontrol_march_conifer_correct_priors.txt – need to download BAMM but then use this to run it
lesfinally2018.tre – Newick formatted tree from Leslie et al. 2018, need to run BAMM

Bamm Script for Publication.R – post divcontrol text run, use this to see results.

Outputs in publication (use Bamm Script for Publication.R)

event_data_man6-8.txt – output file needed to run R script for results in publication
mcmc_out_man6-8.txt – other output file needed to run R script for results in publication

HiSSE

Hisse ultra Final Code 08202024.R – use this to run program
conif_species_fire_trait_use_for_hisse_post_res_redo_20_aug_2024.csv - trait file for hisse
lesfinally2018.tre – Newick formatted tree from Leslie et al. 2018

MiSSE

FireMiSSE.R – use this to run everything
lesfinally2018.tre – Newick formatted tree from Leslie et al. 2018

LSBDS - need revbayes

FireLSBDSext.Rev – use this code first to run
SummFireLSBDS_AncStates.Rev – turn LSBDS anc states into MAP tree files
lesfinally2018.tre – Newick formatted tree from Leslie et al. 2018
LSBDS_checks_and_plots.R - R code to plot stuff

Tables:

Table S1: Master File for Citations. If you use any information from this, please cite this publication. First sheet, quick reference for traits. Second sheet, fully written out citations.
Table S2: HiSSE results, all models all results
Table S3: Reported Ancestral State Reconstruction Results for both All Rates Different and Equal Rates Models.
Table S4: Full Pagel Trait Climate Results

Figures:

Ancestral State Recons

Figure S1-4: Ancestral State Reconstruction Results, All Rates Different
Figure S1b-4b: Ancestral State Reconstruction Results, Equal Rates
Figure S5-8: Jin et al 2021 Grass Stage Ancestral State Reconstruction Plots. 4 plots for both 4 and 16 fossil calibrations and also for the ARD and ER models.
Figure S6-8: Zoomed in BAMM plots

Pagel Climate Trait Relationships

Figure S9: Pagel Climate - Trait side by side phylogenetic plots

BAMM

Note runs 6-8 refer to our final three mcmc runs whose data are in the publication. To keep naming conventions straight, we are keeping these run names

Figure S10-12:Best Clade Shift plots (runs 6,7,8)
Figure S13-15: Zoomed in outpout net diversification plots (runs 6,7,8)
Figure S16-18: Credible Clade Shift Plots (runs 6,7,8)
Figure S19-21: Maximum Shift Credibility plots (runs 6,7,8)

MiSSE

Figure S22: MiSSE Output

LSBDS

Figure S23: LSBDS net diversification plot (effectively a zoomed in version of Figure 3)
Figure S24: Number of shifts plot
Figure S25: LSBDS Lambda
Figure S26: LSBSS mu