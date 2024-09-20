require(phytools)
require(ape)
require(TreeSim)
require(geiger)

##Load Data

getwd()

setwd("/Users/danturck/Documents/Systematics Conifer/Working/")

#get excel file
X<-read.csv("conif_species_fire_trait_post_res_redo_20_aug_2024.csv",row.names=1)

       conif.tree<-read.tree("lesfinally2018.tre")
       conif.tree

       #conif.tree <- ladderize(conif.tree)
       #extraneous in this code: eliminate polytomies (bi.conif meaning binary conifer tree)
       bi.conif<-conif.tree


       plotTree(bi.conif,fsize=0.1,ftype="i",lwd=0.5)

       #make a root
       #bbi.conif<-bind.tip(bi.conif, "rootTip", edge.length=0,where=579)
       #plotTree(bbi.conif, fsize=0.1, lwd=1)

       bi.conif<-multi2di(bi.conif)
       
      

#this sets the traits. In line 17, X[,_] changes the traits 1= bark, 2= serotiny, #3= grass, 4= resprouting

nyDensityMap <- list()
myDensityMap <- list()

for (i in 1:4) {
       bark.adapt <- c()
       bark.adapt<-setNames(X[, i],rownames(X))
       bark.adapt

#make root state fixed 
       #y<-"yes"
       #names(y)<-"rootTip"
       #bbark.adapt<-c(as.character(bark.adapt), y)
       #bbark.adapt<-as.factor(bbark.adapt)
       #names(bbark.adapt)<-c(names(bark.adapt), names(y))
       #bbark.adapt

#give traits distinct colors
       cols<-setNames(c("blue", "green"),levels(bark.adapt))

       levels(bark.adapt)
       bark.adapt<-bark.adapt


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

       

       nodelabels(pie=fit3$lik.anc, cex=0.1)
       #tiplabels(pie=to.matrix(bark.adapt,seq=letters[1:3]),cex=0.5)

       #poop

       nodelabels(pie=fit4$lik.anc, cex=0.1)
       #tiplabels(pie=to.matrix(bark.adapt,seq=letters[1:3]),cex=0.5)

       #chk<-name.check(bi.conif,bark.adapt) 

       #Make simmaps
       #BELOW HERE TREES REBEL

       vtree<-make.simmap(bi.conif,bark.adapt,model="ER", pi='estimated', lwd=1)

       vtrees<-make.simmap(bi.conif,bark.adapt,model="ER",nsim=100)


       mtree<-make.simmap(bi.conif,bark.adapt,model="ARD", pi='estimated', lwd=1)

       mtrees<-make.simmap(bi.conif,bark.adapt,model="ARD",nsim=100)

       #trees<-make.simmap(new.tree,new.data,nsim=100)

nyDensityMap[[i]]<-densityMap(vtrees,states=levels(bark.adapt)[1:2],plot=FALSE)

myDensityMap[[i]]<-densityMap(mtrees,states=levels(bark.adapt)[1:2],plot=FALSE)

       #Summary of ER Simmap
       nbj<-summary(vtrees)
       nbj

       #Summary of ARD Simmap
       mbj<-summary(mtrees)
       mbj

}

par(mfrow=c(2,2))


for (i in 1:4) {
       #Plot ER tree
       plot(nyDensityMap[[i]],fsize=c(0.1,0.1), lwd=1, fsize=1)
       obj<-geo.legend()
}

windows()
par(mfrow=c(2,2))

for (i in 1:4) {

       #Plot ARD tree
       plot(myDensityMap[[i]],colors=cols,fsize=c(0.1,0.1), lwd=1)
       obj<-geo.legend()
}

