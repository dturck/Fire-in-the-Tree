library(hisse)
library(phytools)

#read in data
getwd()

setwd("/Users/danturck/Documents/Systematics Conifer/Working/Hisse_post_ajb/")

con.data<-read.csv("conif_species_fire_trait_use_for_hisse_post_res_redo_20_aug_2024.csv",row.names=1)

#set your trait via column names in csv file
hd<-data.frame(taxa=rownames(con.data), pyro=con.data[,"More"])
head(hd)

gt<-read.tree("lesfinally2018.tre")
gt<-force.ultrametric(gt)


#4x4 rates matrix for hisse no constraints

rates.4x4.free<-TransMatMakerHiSSE(hidden.traits=1) 
rates.4x4.free

#4x4 rates matrix = Null

rates.4x4.equal<-TransMatMakerHiSSE(hidden.traits=1, make.null = TRUE)
rates.4x4.equal

#4x4 rates matrix = 1 

rates.4x4.1<-rates.4x4.free 
rates.4x4.1[!is.na(rates.4x4.1)]<-1 
rates.4x4.1

#8x8 rates no constraints

rates.8x8.free<-TransMatMakerHiSSE(hidden.traits=3) 
rates.8x8.free

#8x8 rates matrix = Null

rates.8x8.equal<-TransMatMakerHiSSE(hidden.traits=3, make.null = TRUE)
rates.8x8.equal

#4x4 rates matrix = 1 

rates.8x8.1<-rates.8x8.free 
rates.8x8.1[!is.na(rates.8x8.1)]<-1 
rates.8x8.1

#define binary rates matrix

rates.bisse<-TransMatMakerHiSSE(hidden.traits=0) 
rates.bisse

#start a binary run

bisse.hmle<-hisse(gt,hd,turnover=c(1,2), eps=c(1,2),hidden.states=FALSE, 
                  trans.rate=rates.bisse)

#inspect bisse model

bisse.hmle

#cid.mle

cid.mle<-hisse(gt,hd,turnover=c(1,1), eps=c(1,1),
               hidden.states=FALSE, trans.rate=rates.bisse)

#new magic function

repar.bd<-function(object,k=2){ 
  pars<-object$solution 
  tt<-pars[grep("turnover",names(pars))][1:k] 
  ee<-pars[grep("eps",names(pars))][1:k] 
  lambda<-tt/(1+ee)
  mu<-tt-lambda 
  nn<-sapply(strsplit(names(tt),"turnover"),function(x) x[2])
  matrix(c(lambda,mu),k,2,dimnames=list(nn,c("lambda","mu")))
}

repar.bd(bisse.hmle)

#single character independent model

cid2.mle<-hisse(gt,hd,f=c(1,1), 
                turnover=c(1,1,2,2),eps=c(1,1,2,2), 
                hidden.states=TRUE,trans.rate=rates.4x4.1)

#Hisse unconstrained

hisse.mle<-hisse(gt,hd,f=c(1,1),
                 hidden.states=TRUE,
                 turnover=c(1,2,3,4),
                 eps=c(1,2,3,4),
                 trans.rate=rates.4x4.free)

#Hisse equal extinction

hisse.ee<-hisse(gt,hd,f=c(1,1),
                hidden.states=TRUE,
                turnover=c(1,2,3,4),
                eps=c(1,1,1,1),
                trans.rate=rates.4x4.free)


#Hisse transitions equal, with Null

hisse.te.null<-hisse(gt,hd,f=c(1,1),
                     hidden.states=TRUE,
                     turnover=c(1,2,3,4),
                     eps=c(1,2,3,4),
                     trans.rate=rates.4x4.equal)

#Hisse transitions equal, with 1x1

hisse.te.1x1<-hisse(gt,hd,f=c(1,1),
                    hidden.states=TRUE,
                    turnover=c(1,2,3,4),
                    eps=c(1,2,3,4),
                    trans.rate=rates.4x4.1)

#Hisse all rates equal, with null

hisse.ae.null<-hisse(gt,hd,f=c(1,1),
                     hidden.states=TRUE,
                     turnover=c(1,2,3,4),
                     eps=c(1,1,1,1),
                     trans.rate=rates.4x4.equal)

#Hisse all rates equal, with 1x1

hisse.ae.1x1<-hisse(gt,hd,f=c(1,1),
                    hidden.states=TRUE,
                    turnover=c(1,2,3,4),
                    eps=c(1,1,1,1),
                    trans.rate=rates.4x4.1)

#CID4 no constraints

cid4.mle<-hisse(gt,hd,f=c(1,1),
                turnover=c(1,1,2,2,3,3,4,4),
                eps=c(1,1,2,2,3,3,4,4),
                hidden.states=TRUE,
                trans.rate=rates.8x8.free)


#CID4 extinction rates equal

cidee.mle<-hisse(gt,hd,f=c(1,1),
                 turnover=c(1,1,2,2,3,3,4,4),
                 eps=c(1,1,1,1,1,1,1,1),
                 hidden.states=TRUE,
                 trans.rate=rates.8x8.free)

#CID4 transition rates equal, with null

cidte.mle.null<-hisse(gt,hd,f=c(1,1),
                      turnover=c(1,1,2,2,3,3,4,4),
                      eps=c(1,1,2,2,3,3,4,4),
                      hidden.states=TRUE,
                      trans.rate=rates.8x8.equal)

#CID4 transition rates equal, with 1x1

cidte.mle.1x1<-hisse(gt,hd,f=c(1,1),
                     turnover=c(1,1,2,2,3,3,4,4),
                     eps=c(1,1,2,2,3,3,4,4),
                     hidden.states=TRUE,
                     trans.rate=rates.8x8.1)


#CID4 all rates equal, with Null

cidae.mle.null<-hisse(gt,hd,f=c(1,1),
                      turnover=c(1,1,2,2,3,3,4,4),
                      eps=c(1,1,1,1,1,1,1,1),
                      hidden.states=TRUE,
                      trans.rate=rates.8x8.equal)

#CID4 all rates equal, with 1x1

cidae.mle.1x1<-hisse(gt,hd,f=c(1,1),
                     turnover=c(1,1,2,2,3,3,4,4),
                     eps=c(1,1,1,1,1,1,1,1),
                     hidden.states=TRUE,
                     trans.rate=rates.8x8.1)


#make function for AIC stuff

logLik.hisse.fit<-function(x,...){
  lik<-x$loglik
  attr(lik,"df")<-(x$AIC+2*lik)/2
  lik
}

#This may fail the first time if running as a block of code. If so just rerun this again.

#make giant dataframe

bark.dataframe<-data.frame(
  model=c("cid, Character independent", "BiSSE", "Char independent with 2 rates","HiSSE with BiSSE, all params free", "HiSSE with BiSSE, extinctions equal", "HiSSE with BiSSE, transitions equal with null", "HiSSE with BiSSE, transitions equal with 1x1", "HiSSE with BiSSE, everything equal with null", "HiSSE with BiSSE, everything equal with 1x1", "Char independent with 4 rates, all params free", "Char independent with 4 rates, extinctions equal", "Char independent with 4 rates, transitions equal with null", "Char independent with 4 rates, transitions equal with 1x1","Char independent with 4 rates, everything equal with null", "Char independent with 4 rates, everything equal with 1x1"),
  
  logL=sapply(list(cid.mle, bisse.hmle, cid2.mle, hisse.mle, hisse.ee, hisse.te.null, hisse.te.1x1, hisse.ae.null, hisse.ae.1x1, cid4.mle, cidee.mle, cidte.mle.null, cidte.mle.1x1, cidae.mle.null, cidae.mle.1x1),logLik),
  
  
  k=sapply(list(cid.mle, bisse.hmle, cid2.mle, hisse.mle, hisse.ee, hisse.te.null, hisse.te.1x1, hisse.ae.null, hisse.ae.1x1, cid4.mle, cidee.mle, cidte.mle.null, cidte.mle.1x1, cidae.mle.null, cidae.mle.1x1),function(x) attr(logLik(x),"df")),
  
  AIC=aic<-sapply(list(cid.mle, bisse.hmle, cid2.mle, hisse.mle, hisse.ee, hisse.te.null, hisse.te.1x1, hisse.ae.null, hisse.ae.1x1, cid4.mle, cidee.mle, cidte.mle.null, cidte.mle.1x1, cidae.mle.null, cidae.mle.1x1),
                  AIC),
  Akaike.weight=unclass(aic.w(aic))
)

#make master table

cid.sol<-(t(cid.mle$solution))
bisse.sol<-(t(bisse.hmle$solution))
cid2.sol<-(t(cid2.mle$solution))
hisse.sol<-(t(hisse.mle$solution)) 
hisse.ee.sol<-(t(hisse.ee$solution))
hisse.te.null.sol<-(t(hisse.te.null$solutio)) 
hisse.te.1x1.sol<-(t(hisse.te.1x1$solution))
hisse.ae.null.sol<-(t(hisse.ae.null$solution))
hisse.ae.1x1.sol<-(t(hisse.ae.1x1$solution))
cid4.sol<-(t(cid4.mle$solution))
cidee.sol<-(t(cidee.mle$solution))
cidte.null.sol<-(t(cidte.mle.null$solution)) 
cidte.1x1.sol<-(t(cidte.mle.1x1$solution))
cidae.null.sol<-(t(cidae.mle.null$solution))
cidae.1x1.sol<-(t(cidae.mle.1x1$solution))

master.rates.table <- rbind(cid.sol, bisse.sol, cid2.sol, hisse.sol, hisse.ee.sol, hisse.te.null.sol, hisse.te.1x1.sol, hisse.ae.null.sol, hisse.ae.1x1.sol, cid4.sol, cidee.sol, cidte.null.sol, cidte.1x1.sol, cidae.null.sol, cidae.1x1.sol)

newcols<-c("turnover0A", "turnover1A", "turnover0B", "turnover1B", "turnover0C", "turnover1C", "turnover0D", "turnover1D", "eps0A", "eps1A", "eps0B", "eps1B", "eps0C", "eps1C", "eps0D", "eps1D", "q0A1A", "q1A0A", "q0B1B", "q1B0B", "q0A0B", "q0B0A", "q1A1B", "q1B1A", "q0C1C", "q1C0C", "q0D1D", "q1D0D", "q0C0D", "q0D0C", "q1C1D", "q1D1C", "q0A0C", "q0C0A", "q1A1C", "q1C1A", "q0B0C", "q0C0B", "q1B1C", "q1C1B", "q0A0D", "q0D0A", "q1A1D", "q1D1A", "q0B0D", "q0D0B", "q1B1D", "q1D1B")

master.rates.table.new<-(master.rates.table[, newcols])

write.csv(master.rates.table.new, "more_Aug_20_2024.csv")

write.csv(bark.dataframe, "more_dataframe_20082024.csv")

echo("poop")
