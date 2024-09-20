require(phytools)
getwd()

#set your directory
setwd("")
conif.tree <- read.tree("lesfinally2018.tre")

#plotTree(conif.tree, fsize = 0.5)


csv.file <- read.csv("climate_csv_07112022_postAJB_res_more_USE.csv",row.names=1,header=TRUE)



#do this first
oeal.list<- list()
#for each climate
for(j in c(1:14)){
  climate<-setNames(csv.file[,j],rownames(csv.file))
  #for each trait do
  for(k in c(23:24)){
    trait<-setNames(csv.file[,k],rownames(csv.file))
    
    fit.pagel<- fitPagel(conif.tree,climate,trait)
    
    clim.name<- colnames(csv.file[j])
    trait.name <- colnames(csv.file[k])
    
    List.name <- paste0(clim.name, " ", trait.name)
    
    oeal.list[[List.name]] <- fit.pagel
    
    print(paste0("done with",clim.name,trait.name))
    
    
  }
  
  
  
}

#Now you will have oeal.list. The insane loop we made for actually extracting data from this is lost.....

#BUT, to look up results from each climate and trait you just call oeal.list[[number in list of climate and trait]]

#The numbered list is below. So if you want to look up say Tropical Savannah Climate and Pyrophilic Trait you would call:

oeal.list[[11]]

#if you want to run each climate and trait individually see Pagel_Climate_Fire_Trait.R in this supplement

#1	Tropical_RF Bark
#2	Tropical_RF Ser
#3	Tropical_RF Grass
#4	Tropical_RF Resprout
#5	Tropical_RF Pyrophillic
#6	Tropical_RF More_than_one_adaptation
#7	Tropical_Monsoon Bark
#8	Tropical_Monsoon Ser
#9	Tropical_Monsoon Grass
#10	Tropical_Monsoon Resprout
#11	Tropical_Monsoon Pyrophillic
#12	Tropical_Monsoon More_than_one_adaptation
#13	Tropical_Savannah Bark
#14	Tropical_Savannah Ser
#15	Tropical_Savannah Grass
#16	Tropical_Savannah Resprout
#17	Tropical_Savannah Pyrophillic
#18	Tropical_Savannah More_than_one_adaptation
#19	Humid_Subtropical Bark
#20	Humid_Subtropical Ser
#21	Humid_Subtropical Grass
#22	Humid_Subtropical Resprout
#23	Humid_Subtropical Pyrophillic
#24	Humid_Subtropical More_than_one_adaptation
#25	Humid_sub_Mont Bark
#26	Humid_sub_Mont Ser
#27	Humid_sub_Mont Grass
#28	Humid_sub_Mont Resprout
#29	Humid_sub_Mont Pyrophillic
#30	Humid_sub_Mont More_than_one_adaptation
#31	Desert Bark
#32	Desert Ser
#33	Desert Grass
#34	Desert Resprout
#35	Desert Pyrophillic
#36	Desert More_than_one_adaptation
#37	Steppe Bark
#38	Steppe Ser
#39	Steppe Grass
#40	Steppe Resprout
#41	Steppe Pyrophillic
#42	Steppe More_tha n_one_adaptation
#43	Mediterranean Bark
#44	Mediterranean Ser
#45	Mediterranean Grass
#46	Mediterranean Resprout
#47	Mediterranean Pyrophillic
#48	Mediterranean More_than_one_adaptation
#49	Med_Cool_sum Bark
#50	Med_Cool_sum Ser
#51	Med_Cool_sum Grass
#52	Med_Cool_sum Resprout
#53	Med_Cool_sum Pyrophillic
#54	Med_Cool_sum More_than_one_adaptation
#55	Oceanic Bark
#56	Oceanic Ser
#57	Oceanic Grass
#58	Oceanic Resprout
#59	Oceanic Pyrophillic
#60	Oceanic More_than_one_adaptation
#61	Cont_Mont Bark
#62	Cont_Mont Ser
#63	Cont_Mont Grass
#64	Cont_Mont Resprout
#65	Cont_Mont Pyrophillic
#66	Cont_Mont More_than_one_adaptation
#67	Humid_Cont_Boreal_subpolar Bark
#68	Humid_Cont_Boreal_subpolar Ser
#69	Humid_Cont_Boreal_subpolar Grass
#70	Humid_Cont_Boreal_subpolar Resprout
#71	Humid_Cont_Boreal_subpolar Pyrophillic
#72	Humid_Cont_Boreal_subpolar More_than_one_adaptation
#73	Mediterranean_Combo Bark
#74	Mediterranean_Combo Ser
#75	Mediterranean_Combo Grass
#76	Mediterranean_Combo Resprout
#77	Mediterranean_Combo Pyrophillic
#78	Mediterranean_Combo More_than_one_adaptation
#79	Cool_Continental Bark
#80	Cool_Continental Ser
#81	Cool_Continental Grass
#82	Cool_Continental Resprout
#83	Cool_Continental Pyrophillic
#84	Cool_Continental More_than_one_adaptation

 
