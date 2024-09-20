library(phytools)

#read in csv

X<-read.csv("climate_csv_07112022_postAJB_res_more_USE.csv",row.names=1,header=TRUE)

bark.adapt<-setNames(X[,1],rownames(X))
bark.adapt

#read in tree
conif.tree<-read.tree("lesfinally2018.tre")
conif.tree

#set climate (X$Climate name in csv), example below
x<-setNames(X$Mediterranean_all_summer_variations,rownames(X))

#set trait (X$Trait name in csv), example below
z<-setNames(X$Resprout,rownames(X))

#get pagel correlation
fit.xz<-fitPagel(conif.tree,x,z)

#get results!
fit.xz

#plot if you wish
plot(fit.xz)
fit.xy