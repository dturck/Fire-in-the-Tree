require(phytools)
require(ape)
require(TreeSim)
require(geiger)

##Load Data

setwd("~/Documents/Systematics Conifer/Working/")

#get excel file
#trait file
X<-read.csv("Pinus_fire_for_ace_with_outgroups.csv",row.names=1)
#X<-read.csv("conif_species_fire_trait_post_res_redo_20_aug_2024.csv",row.names=1)



#conif.tree<-read.tree("lesfinally2018.tre")

#Jin et al has two time calibrated phylogenies. One with 4 calibrations one with 16
conif.tree<-read.tree("cds_time_4_fossils_newick.tre")

#conif.tree<-drop.tip(conif.tree,tip)

#write.tree(conif.tree,"Jin_et_al_2021_Pinus_no_outgroups.tre")

plotTree(conif.tree)

#tip<-c("arg_S001","bre_S003","abi_S002","smi_S005")


##Load Data



#this sets the traits. In line 37, X[,_] changes the traits 4= bark, 5= serotiny
#6= grass, 7= resprouting 8= pyrophillic
bark.adapt<-setNames(X[,6],rownames(X))
bark.adapt


#extraneous in this code: eliminate polytomies (bi.conif meaning binary conifer tree)
bi.conif<-conif.tree


plotTree(bi.conif,fsize=0.1,ftype="i",lwd=1)

#make a root
#bbi.conif<-bind.tip(bi.conif, "rootTip", edge.length=0,where=579)
#plotTree(bbi.conif, fsize=0.1, lwd=1)
#
#bi.conif<-multi2di(bbi.conif)

#bi.conif<-bbi.conif 

#make root state fixed at no
#y<-"no"
#names(y)<-"rootTip"
#bbark.adapt<-c(as.character(bark.adapt), y)
#bbark.adapt<-as.factor(bbark.adapt)
#names(bbark.adapt)<-c(names(bark.adapt), names(y))
#bbark.adapt

#give traits distinct colors
cols<-setNames(c("blue", "green"),levels(bark.adapt))

levels(bark.adapt)
bark.adapt<-bbark.adapt


#place char state on the root tip and current char states on all other tips
tiplabels(pie=to.matrix(bark.adapt[bi.conif$tip.label],
                        levels(bark.adapt)),piecol=cols,cex=0.1)

#ancestral char state equalrates
dsb<-bi.conif
dsb$edge.length[dsb$edge.length==0]<-max(nodeHeights(bi.conif))*1e-6
fit3<-ace(bark.adapt,dsb,type="discrete",model="ER")
fit3

#aic score for ER
aic3<-AIC(fit3, k=2)
aic3  


#ancestral char states all rates different
dsb<-bi.conif
dsb$edge.length[dsb$edge.length==0]<-max(nodeHeights(bi.conif))*1e-6
fit4<-ace(bark.adapt,dsb,type="discrete",model="ARD")
fit4

#aic score for ARD
aic4<-AIC(fit4, k=2)
aic4

#dsb<-bi.conif
#dsb$edge.length[dsb$edge.length==0]<-max(nodeHeights(bi.conif#))*1e-6
#fit5<-ace(bark.adapt,dsb,type="discrete",model="SYM")
#fit5
#
##aic score for ARD
#aic5<-AIC(fit5, k=2)
#aic5
#
#dsb<-bi.conif
#dsb$edge.length[dsb$edge.length==0]<-max(nodeHeights(bi.conif#))*1e-6
#fit6<-ace(bark.adapt,dsb,type="discrete",model="IR")
#fit6
#
##aic score for ARD
#aic6<-AIC(fit6, k=2)
#aic6
#
#nodelabels(pie=fit3$lik.anc, cex=0.1)
##tiplabels(pie=to.matrix(bark.adapt,seq=letters[1:3]),cex=0.5)
#
#poop
#
#nodelabels(pie=fit4$lik.anc, cex=0.1)
#tiplabels(pie=to.matrix(bark.adapt,seq=letters[1:3]),cex=0.5)




vtree<-make.simmap(bi.conif,bark.adapt,model="ER", pi='estimated', lwd=1)

vtrees<-make.simmap(bi.conif,bark.adapt,model="ER",nsim=100)


mtree<-make.simmap(bi.conif,bark.adapt,model="ARD", pi='estimated', lwd=1)

mtrees<-make.simmap(bi.conif,bark.adapt,model="ARD",nsim=100)

#trees<-make.simmap(new.tree,new.data,nsim=100)

nyDensityMap<-densityMap(vtrees,states=levels(bark.adapt)[1:2],plot=FALSE)
myDensityMap<-densityMap(mtrees,states=levels(bark.adapt)[1:2],plot=FALSE, colors=cols)

#Plot ER tree
plot(nyDensityMap,fsize=c(0.1,0.1), lwd=1, fsize=1)
obj<-geo.legend()

#Plot ARD tree
plot(myDensityMap,colors=cols,fsize=c(0.1,0.1), lwd=1)
obj<-geo.legend()


#Summary of ER Simmam
nbj<-summary(vtrees)
nbj

#Summary of ARD Simmam
mbj<-summary(mtrees)
mbj



