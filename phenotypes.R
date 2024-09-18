setwd("C:/Users/Steph/OneDrive - Texas A&M University/Salmonids/data")

###Read in packages

pkgs<-c("ggplot2","metafor","metaDigitise","tidyverse","gridExtra","matrixcalc",
        "formattable","ape","raster","sp","MuMIn","ggpubr",
        "glmulti","brms")

lapply(pkgs,library,character.only=T); rm(pkgs)
source("C:/Users/Steph/OneDrive - Texas A&M University/Salmonids/salmonids/functions_phenotypes.R")

theme_set(theme_classic()) #set theme for graphing

###Colours
col.sal<-"#B4674E"
col.cor<-"#2E8289"

###Read in data

df<-read.csv("lakes_summary.csv",stringsAsFactors=F)
m<-read.csv("morphs_summary.csv",stringsAsFactors=F)
morphs<-read.csv("morphs.csv",stringsAsFactors=F,fileEncoding="latin1")

morphs<-subset(morphs,morphs$Species!="Prosopium_coulterii")
m<-subset(m,m$Species!="Prosopium_coulterii")
df<-subset(df,df$Lake%in%m$Lake)


df$X<-NULL; morphs$X<-NULL

#get log variable for lake size to normalize (somewhat)
df$logArea<-log(df$Area)
df$logDepthMax<-log(df$Depth_max)



#get climatic variables from worldClim for temperature and precipitation
r <- getData("worldclim",var="bio",res=10)[[c(1,4,12,15)]]

#Temp = Temperature, Prec = Precipitation, _sea = seasonality, i.e. coef of var
names(r) <- c("Temp_mean","Temp_sea","Prec_mean","Prec_sea")

#extract climatic vars for each lake based on latitude and longitude, attach to dataframe
df<-cbind(df,raster::extract(r,SpatialPoints(cbind(x=df$Lon,y=df$Lat),proj4string=r@crs)))
hist(log(df$Prec_sea),breaks=30)
df$Temp_sea<-df$Temp_sea/100
df.pc<-prcomp(df[,19:22])

#add squared values to test if relationships non-linear
df$logArea.2<-df$logArea^2; df$logDepthMax.2<-df$logDepthMax^2; df$Lat.2<-df$Lat^2
df$Prec.2<-df$Prec_mean^2; df$Temp.2<-df$Temp_mean^2

df$Continent<-df$Country
df$Continent[df$Continent%in%c("US","Canada","Greenland")]<-"N_America"
df$Continent[df$Continent%in%c("Iceland","UK","Switzerland","Norway","Finland",
                               "Finland/Sweden","Germany")]<-"Europe"
df$Continent[df$Continent=="Russia"]<-"Asia"
unique(df$Continent) #check levels are N America, Europe, and Asia

#relative depth - higher for deep lakes with small surface area
#use on log scale
df$Zr<-50*df$Depth_max*(sqrt(pi)/sqrt(df$Area))

df1<-df%>%mutate(logZr=log(Zr))%>%
  dplyr::select(Lake,logArea,logZr,logDepthMax,Temp_mean,Temp_sea,Prec_mean,Prec_sea,Lat)%>%
  drop_na()
corMat<-cor(df1%>%dplyr::select(-Lake))%>%as.data.frame()%>%
  rownames_to_column(var="lakeVar1")%>%
  pivot_longer(-lakeVar1,names_to="lakeVar2",values_to="corCoef")

ggplot(corMat,aes(x=lakeVar1,y=lakeVar2,fill=corCoef))+
  geom_tile()+coord_equal()+
  xlab("Lake feature")+ylab("Lake feature")+
  scale_fill_gradient2(name="Correlation\ncoefficient")

df1.pca<-prcomp(df1%>%dplyr::select(-Lake),scale=T)
df1.pc<-df1.pca$rotation%>%as.data.frame()%>%
  rownames_to_column(var="lakeVar")%>%
  pivot_longer(-lakeVar,names_to="PC",values_to="loading")
summary(df1.pca)

#PC1 = climate, PC2 = lake size, PC3 = lake ratio, PC4 = latitude
ggplot(df1.pc%>%filter(PC%in%paste("PC",1:4,sep="")),
       aes(x=PC,y=lakeVar,fill=loading))+
  geom_tile()+coord_equal()+
  scale_fill_gradient2(high="darkslategrey",low="lightpink4")+
  xlab("PC axis")

df1.pca.x<-df1.pca$x%>%as.data.frame()%>%
  dplyr::select(PC1,PC2,PC3,PC4)%>%
  rename("PC1_climate"=PC1,"PC2_size"=PC2,
         "PC3_depthRatio"=PC3,"PC4_Latitude"=PC4)%>%
  mutate(Lake=df1$Lake)

df<-df%>%left_join(df1.pca.x)

#directory for figures output
setwd("C:/Users/Steph/OneDrive - Texas A&M University/Salmonids/figures/phenotypes")

###Get effect sizes

#fix problems with zero error by removing measure from analysis
morphs$GRc_e[morphs$GRc_e==0]<-NA

traitColNames<-c("GRc_mean","GRc_n","GRc_e","GRc_e_type",
                 "length_mean","length_n","length_e","length_e_type","length_units","length_type",
                 "age","age_n","age_e","age_e_type",
                 "lengthAsymp_mean","lengthAsymp_n","lengthAsymp_e","lengthAsymp_e_type","lengthAsymp_type",
                 "N15_mean","N15_n","N15_e","N15_e_type",
                 "C13_mean","C13_n","C13_e","C13_e_type",
                 "UpperJaw_mean","UpperJaw_n","UpperJaw_e","UpperJaw_e_type",
                 "LowerJaw_mean","LowerJaw_n","LowerJaw_e","LowerJaw_e_type",
                 "Gape_mean","Gape_n","Gape_e","Gape_e_type",
                 "Eye_mean","Eye_n","Eye_e","Eye_e_type",
                 "HeadLength_mean","HeadLength_n","HeadLength_e","HeadLength_e_type",
                 "Snout_mean","Snout_n","Snout_e","Snout_e_type",
                 "PectoralFin_mean","PectoralFin_n","PectoralFin_e","PectoralFin_e_type",
                 "DorsalFin_mean","DorsalFin_n","DorsalFin_e","DorsalFin_e_type",
                 "CaudalFin_mean","CaudalFin_n","CaudalFin_e","CaudalFin_e_type",
                 "AnalFin_mean","AnalFin_n","AnalFin_e","AnalFin_e_type",
                 "PelvicFin_mean","PelvicFin_n","PelvicFin_e","PelvicFin_e_type",
                 "Depth","Depth_n","Depth_e","Depth_e_type",
                 "BodyPC1_mean","BodyPC1_n","BodyPC1_e","BodyPC1_e_type",
                 "HeadPC1_mean","HeadPC1_n","HeadPC1_e","HeadPC1_e_type")

colNums=which(colnames(morphs)%in%c("GRc_mean","GRc_n","GRc_e","GRc_e_type"))
#use counts for gill rakers, not lengths, because more data
dat.gr<-calc.dat.tr(colNums,"gr_count") #gill rakers adjacent
dat.gr.o<-calc.dat.tr_outer(colNums,"gr_count") #gill rakers outer
dat.gr.cv<-calc.datCV.tr(colNums,"gr_count") #gill rakers lnCVR

colNums=which(colnames(morphs)%in%c("length_mean","length_n","length_e","length_e_type"))
dat.fl<-calc.dat.tr(colNums,"body_length") #body length adjacent
dat.fl.o<-calc.dat.tr_outer(colNums,"body_length") #body length outer
dat.fl.cv<-calc.datCV.tr(colNums,"body_length") #body length lnCVR

colNums=which(colnames(morphs)%in%c("age","age_n","age_e","age_e_type"))
dat.ag<-calc.dat.tr(colNums,"age")#age adjacent
dat.ag.o<-calc.dat.tr_outer(colNums,"age") #age outer
dat.ag.cv<-calc.datCV.tr(colNums,"age") #age lnCVR

colNums=which(colnames(morphs)%in%c("lengthAsymp_mean","lengthAsymp_n","lengthAsymp_e","lengthAsymp_e_type"))
dat.al<-calc.dat.tr(colNums,"asymptotic_length") #asymptotic length adjacent
dat.al.o<-calc.dat.tr_outer(colNums,"asymptotic_length") #asymptotic length outer
dat.al.cv<-calc.datCV.tr(colNums,"asymptotic_length")

#isotopes
colNums=which(colnames(morphs)%in%c("N15_mean","N15_n","N15_e","N15_e_type"))
dat.N15<-calc.dat.tr(colNums,"N15")
dat.N15.o<-calc.dat.tr_outer(colNums,"N15")
dat.N15.cv<-calc.datCV.tr(colNums,"N15")

colNums=which(colnames(morphs)%in%c("C13_mean","C13_n","C13_e","C13_e_type"))
dat.C13.o<-calc.dat.tr_outer(colNums,"C13")
dat.C13<-calc.dat.tr(colNums,"C13")
dat.C13.cv<-calc.datCV.tr(colNums,"C13")


colNums1=which(colnames(morphs)%in%c("UpperJaw_mean","UpperJaw_n","UpperJaw_e","UpperJaw_e_type"))
colNums2=which(colnames(morphs)%in%c("LowerJaw_mean","LowerJaw_n","LowerJaw_e","LowerJaw_e_type"))
colNums3=which(colnames(morphs)%in%c("Gape_mean","Gape_n","Gape_e","Gape_e_type"))
colNums4=which(colnames(morphs)%in%c("Eye_mean","Eye_n","Eye_e","Eye_e_type"))
colNums5=which(colnames(morphs)%in%c("HeadLength_mean","HeadLength_n","HeadLength_e","HeadLength_e_type"))
colNums6=which(colnames(morphs)%in%c("Snout_mean","Snout_n","Snout_e","Snout_e_type"))
#trophic traits: upper jaw, lower jaw, gape, eye, head length, snout length
#trophic traits adjacent
dat.troph<-rbind(calc.dat.tr(colNums1,"upper_jaw"),calc.dat.tr(colNums2,"lower_jaw"),
                 calc.dat.tr(colNums3,"gape"),calc.dat.tr(colNums4,"eye"),
                 calc.dat.tr(colNums5,"head"),calc.dat.tr(colNums6,"snout"))

#unique(dat.troph$Lake) #only 17 lakes here

#trophic traits outer
dat.troph.o<-rbind(calc.dat.tr_outer(colNums1,"upper_jaw"),calc.dat.tr_outer(colNums2,"lower_jaw"),
                   calc.dat.tr_outer(colNums3,"gape"),calc.dat.tr_outer(colNums4,"eye"),
                   calc.dat.tr_outer(colNums5,"head"),calc.dat.tr_outer(colNums6,"snout"))

dat.troph.cv<-rbind(calc.datCV.tr(colNums1,"upper_jaw"),calc.datCV.tr(colNums2,"lower_jaw"),
                 calc.datCV.tr(colNums3,"gape"),calc.datCV.tr(colNums4,"eye"),
                 calc.datCV.tr(colNums5,"head"),calc.datCV.tr(colNums6,"snout"))

colNums1=which(colnames(morphs)%in%c("PectoralFin_mean","PectoralFin_n","PectoralFin_e","PectoralFin_e_type"))
colNums2=which(colnames(morphs)%in%c("DorsalFin_mean","DorsalFin_n","DorsalFin_e","DorsalFin_e_type"))
colNums3=which(colnames(morphs)%in%c("CaudalFin_mean","CaudalFin_n","CaudalFin_e","CaudalFin_e_type"))
colNums4=which(colnames(morphs)%in%c("AnalFin_mean","AnalFin_n","AnalFin_e","AnalFin_e_type"))
colNums5=which(colnames(morphs)%in%c("PelvicFin_mean","PelvicFin_n","PelvicFin_e","PelvicFin_e_type"))
#fins adjacent
dat.fin<-rbind(calc.dat.tr(colNums1,"pectoral"),calc.dat.tr(colNums2,"dorsal"),
               calc.dat.tr(colNums3,"caudal"),calc.dat.tr(colNums4,"anal"),
               calc.dat.tr(colNums5,"pelvic"))
#fins outer
dat.fin.o<-rbind(calc.dat.tr_outer(colNums1,"pectoral"),calc.dat.tr_outer(colNums2,"dorsal"),
                 calc.dat.tr_outer(colNums3,"caudal"),calc.dat.tr_outer(colNums4,"anal"),
                 calc.dat.tr_outer(colNums5,"pelvic"))

dat.fin.cv<-rbind(calc.datCV.tr(colNums1,"pectoral"),calc.datCV.tr(colNums2,"dorsal"),
               calc.datCV.tr(colNums3,"caudal"),calc.datCV.tr(colNums4,"anal"),
               calc.datCV.tr(colNums5,"pelvic"))


colNums=which(colnames(morphs)%in%c("Depth","Depth_n","Depth_e","Depth_e_type"))
#fix problems with zero error
morphs$Depth_e[morphs$Depth_e==0]<-NA
#depth adjacent 
dat.dep<-calc.dat.tr(colNums,"depth")
#depth outer
dat.dep.o<-calc.dat.tr_outer(colNums,"depth")
dat.dep.cv<-calc.datCV.tr(colNums,"depth")

colNums=which(colnames(morphs)%in%c("BodyPC1_mean","BodyPC1_n","BodyPC1_e","BodyPC1_e_type"))
#PCA body adjacent
dat.pcb<-rbind(calc.dat.tr(colNums,"PC1_body"))#,calc.dat.tr(20:23,"PC2_body"))
#PCA body outer
dat.pcb.o<-rbind(calc.dat.tr_outer(colNums,"PC1_body"))#,calc.dat.tr_outer(20:23,"PC2_body"))
#lnCVR doesn't work with negative and positive values

colNums=which(colnames(morphs)%in%c("HeadPC1_mean","HeadPC1_n","HeadPC1_e","HeadPC1_e_type"))
#PCA head adjacent
dat.pch<-rbind(calc.dat.tr(colNums,"PC1_head"))#,calc.dat.tr(29:32,"PC2_head"))
#PCA head outer
dat.pch.o<-rbind(calc.dat.tr_outer(colNums,"PC1_head"))#,calc.dat.tr_outer(29:32,"PC2_head"))




dat.tr<-rbind(dat.fl,dat.gr,dat.N15,dat.ag,dat.al,dat.C13,
              dat.dep,dat.fin,dat.pcb,dat.pch,dat.troph)
dat.tr.o<-rbind(dat.fl.o,dat.gr.o,dat.N15.o,dat.ag.o,dat.al.o,dat.C13.o,
                dat.dep.o,dat.fin.o,dat.pcb.o,dat.pch.o,dat.troph.o)
dat.tr.cv<-rbind(dat.fl.cv,dat.gr.cv,dat.N15.cv,dat.ag.cv,dat.al.cv,dat.C13.cv,
                dat.dep.cv,dat.fin.cv,dat.troph.cv)



#figure out sample sizes for each trait / genus combo
sampleSizes<-dat.tr%>%group_by(Genus,Trait)%>%summarise(sample_size=n())
ggplot(data=sampleSizes,aes(x=Trait,y=sample_size,fill=Genus))+
  geom_bar(stat='identity')+theme_classic()+
  scale_fill_manual(values=c(col.cor,col.sal))

#find traits that have a large enough sample size for a genus
keepTraits<-sampleSizes%>%filter(sample_size>10)

#only keep traits will large enough sample size
keepTraits$grouping<-paste(keepTraits$Genus,keepTraits$Trait,sep=".")
#in adjacent SMD
dat.tr$grouping<-paste(dat.tr$Genus,dat.tr$Trait,sep=".")
dat.tr<-dat.tr%>%filter(dat.tr$grouping%in%keepTraits$grouping)
dat.tr$grouping<-NULL
#in outer SMD
dat.tr.o$grouping<-paste(dat.tr.o$Genus,dat.tr.o$Trait,sep=".")
dat.tr.o<-dat.tr.o%>%filter(dat.tr.o$grouping%in%keepTraits$grouping)
dat.tr.o$grouping<-NULL
#in lnCVR
dat.tr.cv$grouping<-paste(dat.tr.cv$Genus,dat.tr.cv$Trait,sep=".")
dat.tr.cv<-dat.tr.cv%>%filter(dat.tr.cv$grouping%in%keepTraits$grouping)
dat.tr.cv$grouping<-NULL



#make glmulti function
#use ML because this is a model fitting exercise

rma.glmulti <- function(formula, data, ...)
  rma(formula, vi, data=data, method="ML", ...)


#create sets of predictor variables
#include number of morphs for trait divergence but not 
#abiotic_vars<-c("Lat","logArea","logDepthMax","Prec_mean","Temp_mean","Prec_sea","Temp_sea")
tr_vars<-c("PC1_climate","PC2_size","PC3_depthRatio","PC4_Latitude")

dat.tr.cor<-dat.tr%>%filter(Genus=="Coregonus")%>%
  drop_na(all_of(tr_vars))

dat.tr.sal<-dat.tr%>%filter(Genus=="Salvelinus")%>%
  drop_na(all_of(tr_vars))


#apply function to coregonus body length data
lm.cor.fl<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                data=dat.tr.cor%>%filter(Trait=="body_length"),
                level=1,fitfunction=rma.glmulti,
                crit="aicc",confsetsize=16)

lm.sal.fl<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                data=dat.tr.sal%>%filter(Trait=="body_length"),
                level=1,fitfunction=rma.glmulti,
                crit="aicc",confsetsize=256)

lm.cor.gr<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                data=dat.tr.cor%>%filter(Trait=="gr_count"),
                level=1,fitfunction=rma.glmulti,
                crit="aicc",confsetsize=16)

lm.sal.ag<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                data=dat.tr.sal%>%filter(Trait=="age"),
                level=1,fitfunction=rma.glmulti,
                crit="aicc",confsetsize=256)

#get a list of lakes used in trait analysis
lakes.tr<-c(dat.tr.sal%>%filter(Trait%in%c("age","body_length"))%>%
  pull(Lake)%>%unique(),
dat.tr.cor%>%filter(Trait%in%c("gr_count","body_length"))%>%
  pull(Lake)%>%unique())

#calculate weights for each variable

lm1.tr.weights<-rbind(calc_AICc_weights(lm.cor.fl,"Coregonus","body_length"),
                  calc_AICc_weights(lm.sal.fl,"Salvelinus","body_length"),
                  calc_AICc_weights(lm.cor.gr,"Coregonus","gill_raker_count"),
                  calc_AICc_weights(lm.sal.ag,"Salvelinus","age"))


#refit with strong predictors only
lm1.tr.weights%>%filter(weight_sum>0.7)

lm.sal.fl<-rma(yi,vi,mods=cbind(PC1_climate,PC2_size,PC4_Latitude),
                   data=dat.tr.sal%>%filter(Trait=="body_length"),method="REML")

summary(lm.sal.fl)

lm.cor.fl<-rma(yi,vi,mods=PC4_Latitude,
               data=dat.tr.cor%>%filter(Trait=="body_length"),method="REML")

summary(lm.cor.fl)




p1<-ggplot(data=lm1.tr.weights%>%filter(trait=="body_length"),
           aes(x=pred_var,#x=reorder(pred_var,-weight_sum),
                                y=weight_sum,fill=genus))+
  facet_wrap(~genus,ncol=1)+geom_bar(stat="identity")+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none",text=element_text(size=12))+
  ylim(0,1)+
  xlab("Abiotic factor")+ylab("Relative importance")+
  scale_x_discrete(labels=c("PC1_climate"="PC1\nclimate","PC2_size"="PC2\ndimensions",
                            "PC3_depthRatio"="PC3\ndepth ratio","PC4_Latitude"="PC4\nlatitude"))

dat.tr.fl<-dat.tr%>%filter(Trait=="body_length")
p2<-ggplot(dat.tr.fl,aes(x=as.factor(Lacustrine_morphs),y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor,col.sal))+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("Number of ecotypes")+ylab("Adjacent SMD")+
  coord_fixed(ratio=0.5)

#make figure 4
fig4<-ggarrange(p2,p1,ncol=2,labels=c("A","B"))
fig4
ggsave("Fig4.png",plot=fig4,height=12,width=20,units="cm",bg="white")


s8e<-ggplot(data=lm1.tr.weights%>%filter(trait!="body_length"),
       aes(x=pred_var,
           y=weight_sum,fill=genus))+
  facet_wrap(~genus+trait,ncol=2)+geom_bar(stat="identity")+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none",text=element_text(size=12))+
  ylim(0,1)+
  xlab("Abiotic factor")+ylab("Relative importance")+
  scale_x_discrete(labels=c("PC1_climate"="PC1\nclimate","PC2_size"="PC2\ndimensions",
                            "PC3_depthRatio"="PC3\ndepth ratio","PC4_Latitude"="PC4\nlatitude"))


s8a<-ggplot(dat.tr.cor%>%filter(Trait=="body_length"),
       aes(x=PC4_Latitude,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC4 latitude")+ylab("Adjacent SMD\nbody length")
s8b<-ggplot(dat.tr.sal%>%filter(Trait=="body_length"),
            aes(x=PC1_climate,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.sal))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC1 climate")+ylab("Adjacent SMD\nbody length")
s8c<-ggplot(dat.tr.sal%>%filter(Trait=="body_length"),
            aes(x=PC2_size,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.sal))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC2 dimensions")+ylab("Adjacent SMD\nbody length")
s8d<-ggplot(dat.tr.sal%>%filter(Trait=="body_length"),
            aes(x=PC4_Latitude,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.sal))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC4 latitude")+ylab("Adjacent SMD\nbody length")

figS8<-ggarrange(s8e,ggarrange(s8a,s8b,s8c,s8d),nrow=2,labels=LETTERS,
          heights=c(1.2,2))

ggsave("FigS8.png",plot=figS8,height=20,width=20,units="cm",bg="white")

morph_count<-keepTraits
morph_count$slope.mean<-NA; morph_count$slope.cI<-NA
morph_count$slope.diff<-NA

for(i in 1:nrow(morph_count)){
  lm.count<-lm(yi~Lacustrine_morphs,data=dat.tr%>%
                 filter(Trait==morph_count$Trait[i]&Genus==morph_count$Genus[i]))
  morph_count$slope.cI[i]<-paste(round(confint(lm.count)[2,1],3),
                                 round(confint(lm.count)[2,2],3),sep=", ")
  morph_count$slope.mean[i]<-round(lm.count$coefficients[2],3)
  if(between(0,round(confint(lm.count)[2,1],3),
             round(confint(lm.count)[2,2],3))==F){
    morph_count$slope.diff[i]<-"*" }else{
      morph_count$slope.diff[i]<-""
    }}



##Outer SMD

dat.tr.o.cor<-dat.tr.o%>%filter(Genus=="Coregonus")%>%
  drop_na(all_of(tr_vars))

dat.tr.o.sal<-dat.tr.o%>%filter(Genus=="Salvelinus")%>%
  drop_na(all_of(tr_vars))


#apply function to coregonus body length data
lm.o.cor.fl<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                   data=dat.tr.o.cor%>%filter(Trait=="body_length"),
                   level=1,fitfunction=rma.glmulti,
                   crit="aicc",confsetsize=16)

lm.o.sal.fl<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                   data=dat.tr.o.sal%>%filter(Trait=="body_length"),
                   level=1,fitfunction=rma.glmulti,
                   crit="aicc",confsetsize=256)

lm.o.cor.gr<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                   data=dat.tr.o.cor%>%filter(Trait=="gr_count"),
                   level=1,fitfunction=rma.glmulti,
                   crit="aicc",confsetsize=16)

lm.o.sal.ag<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                   data=dat.tr.o.sal%>%filter(Trait=="age"),
                   level=1,fitfunction=rma.glmulti,
                   crit="aicc",confsetsize=256)


lm1.tr.o.weights<-rbind(calc_AICc_weights(lm.o.cor.fl,"Coregonus","body_length"),
                        calc_AICc_weights(lm.o.sal.fl,"Salvelinus","body_length"),
                        calc_AICc_weights(lm.o.cor.gr,"Coregonus","gill_raker_count"),
                        calc_AICc_weights(lm.o.sal.ag,"Salvelinus","age"))

#refit with strong predictors only
lm1.tr.o.weights%>%filter(weight_sum>0.7)

rma(yi,vi,mods=PC4_Latitude,
    data=dat.tr.o.cor%>%filter(Trait=="body_length"),method="REML")


rma(yi,vi,mods=PC1_climate,
    data=dat.tr.o.sal%>%filter(Trait=="body_length"),method="REML")

rma(yi,vi,mods=PC4_Latitude,
    data=dat.tr.o.cor%>%filter(Trait=="gr_count"),method="REML")

rma(yi,vi,mods=PC1_climate,
    data=dat.tr.o.sal%>%filter(Trait=="age"),method="REML")



s9a<-ggplot(data=lm1.tr.o.weights,
       aes(x=pred_var,
           y=weight_sum,fill=genus))+
  facet_wrap(~genus+trait,ncol=2)+geom_bar(stat="identity")+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none",text=element_text(size=12))+
  ylim(0,1)+
  xlab("Abiotic factor")+ylab("Relative importance")+
  scale_x_discrete(labels=c("PC1_climate"="PC1\nclimate","PC2_size"="PC2\ndimensions",
                            "PC3_depthRatio"="PC3\ndepth ratio","PC4_Latitude"="PC4\nlatitude"))

s9b<-ggplot(dat.tr.o.cor%>%filter(Trait=="body_length"),
       aes(x=PC4_Latitude,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC4 latitude")+ylab("Outer SMD\nbody length")

s9c<-ggplot(dat.tr.o.cor%>%filter(Trait=="gr_count"),
            aes(x=PC4_Latitude,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC4 latitude")+ylab("Outer SMD\ngill rakers")

s9d<-ggplot(dat.tr.o.sal%>%filter(Trait=="body_length"),
            aes(x=PC1_climate,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.sal))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC1 climate")+ylab("Outer SMD\nbody length")

s9e<-ggplot(dat.tr.o.sal%>%filter(Trait=="age"),
            aes(x=PC1_climate,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.sal))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC1 climate")+ylab("Outer SMD\nage")

figS9<-ggarrange(s9a,ggarrange(s9b,s9c,s9d,s9e),nrow=2,labels=LETTERS)

ggsave("FigS9.png",plot=figS9,height=25,width=20,units="cm",bg="white")

morph_count.o<-keepTraits
morph_count.o$slope.mean<-NA; morph_count.o$slope.cI<-NA
morph_count.o$slope.diff<-NA

for(i in 1:nrow(morph_count.o)){
  lm.count<-lm(yi~Lacustrine_morphs,data=dat.tr.o%>%
                 filter(Trait==morph_count.o$Trait[i]&Genus==morph_count.o$Genus[i]))
  morph_count.o$slope.cI[i]<-paste(round(confint(lm.count)[2,1],3),
                                 round(confint(lm.count)[2,2],3),sep=", ")
  morph_count.o$slope.mean[i]<-round(lm.count$coefficients[2],3)
  if(between(0,round(confint(lm.count)[2,1],3),
             round(confint(lm.count)[2,2],3))==F){
    morph_count.o$slope.diff[i]<-"*" }else{
      morph_count.o$slope.diff[i]<-"" }}



##lnCVR

dat.tr.cv.cor<-dat.tr.cv%>%filter(Genus=="Coregonus")%>%
  drop_na(all_of(tr_vars))

dat.tr.cv.sal<-dat.tr.cv%>%filter(Genus=="Salvelinus")%>%
  drop_na(all_of(tr_vars))


#apply function to coregonus body length data
lm.cv.cor.fl<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                     data=dat.tr.cv.cor%>%filter(Trait=="body_length"),
                     level=1,fitfunction=rma.glmulti,
                     crit="aicc",confsetsize=16)

lm.cv.sal.fl<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                     data=dat.tr.cv.sal%>%filter(Trait=="body_length"),
                     level=1,fitfunction=rma.glmulti,
                     crit="aicc",confsetsize=256)

lm.cv.cor.gr<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                     data=dat.tr.cv.cor%>%filter(Trait=="gr_count"),
                     level=1,fitfunction=rma.glmulti,
                     crit="aicc",confsetsize=16)

lm.cv.sal.ag<-glmulti(yi~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
                     data=dat.tr.cv.sal%>%filter(Trait=="age"),
                     level=1,fitfunction=rma.glmulti,
                     crit="aicc",confsetsize=256)

lm1.tr.cv.weights<-rbind(calc_AICc_weights(lm.cv.cor.fl,"Coregonus","body_length"),
                        calc_AICc_weights(lm.cv.sal.fl,"Salvelinus","body_length"),
                        calc_AICc_weights(lm.cv.cor.gr,"Coregonus","gill_raker_count"),
                        calc_AICc_weights(lm.cv.sal.ag,"Salvelinus","age"))


lm1.tr.cv.weights%>%filter(weight_sum>0.7)

rma(yi,vi,mods=PC4_Latitude,
    data=dat.tr.cv.cor%>%filter(Trait=="body_length"),method="REML")

rma(yi,vi,mods=PC3_depthRatio,
    data=dat.tr.cv.sal%>%filter(Trait=="age"),method="REML")


s10a<-ggplot(data=lm1.tr.cv.weights%>%filter(trait=="body_length"),
       aes(x=pred_var,
           y=weight_sum,fill=genus))+
  facet_wrap(~genus+trait,ncol=2)+geom_bar(stat="identity")+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none",text=element_text(size=12))+
  ylim(0,1)+
  xlab("Abiotic factor")+ylab("Relative importance")+
  scale_x_discrete(labels=c("PC1_climate"="PC1\nclimate","PC2_size"="PC2\ndimensions",
                            "PC3_depthRatio"="PC3\ndepth ratio",
                            "PC4_Latitude"="PC4\nlatitude"))

s10b<-ggplot(dat.tr.o.sal%>%filter(Trait=="body_length"),
             aes(x=PC4_Latitude,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC4_latitude")+ylab("lnCVR\nbody length")

s10c<-ggplot(dat.tr.o.sal%>%filter(Trait=="age"),
            aes(x=PC3_depthRatio,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.sal))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC3 depth ratio")+ylab("lnCVR\nage")

figS10<-ggarrange(s10a,ggarrange(s10b,s10c),nrow=2,heights=c(1.2,1),labels=LETTERS)

ggsave("FigS10.png",plot=figS10,height=16,width=20,units="cm",bg="white")

lm1.tr.weights$estimate<-"adjacent SMD"
lm1.tr.o.weights$estimate<-"outer SMD"
lm1.tr.cv.weights$estimate<-"lnCVR"

lm1.weights<-rbind(lm1.tr.weights,lm1.tr.o.weights,lm1.tr.cv.weights)
lm1.weights<-lm1.weights%>%arrange(estimate,trait,genus,pred_var)
#write.csv(lm1.weights,"trait_weights.240716.csv")


p1<-ggplot(data=lm1.tr.cv.weights%>%filter(trait!="body_length"),
       aes(x=pred_var,
           y=weight_sum,fill=genus))+
  facet_wrap(~genus+trait,ncol=1)+geom_bar(stat="identity")+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none",text=element_text(size=12))+
  ylim(0,1)+
  xlab("Abiotic factor")+ylab("Relative importance")+
  scale_x_discrete(labels=c("PC1_climate"="PC1\nclimate","PC2_size"="PC2\ndimensions",
                            "PC3_depthRatio"="PC3\ndepth ratio",
                            "PC4_Latitude"="PC4\nlatitude"))

p2<-ggplot(dat.tr.cv.cor%>%filter(Trait=="gr_count"),
       aes(x=PC1_climate,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC1 - climate")+ylab("lnCVR - gill rakers")

p3<-ggplot(dat.tr.cv.sal%>%filter(Trait=="age"),
       aes(x=PC3_depthRatio,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.sal))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("PC3 - depth ratio")+ylab("lnCVR - age")

fig5<-ggarrange(ggarrange(p2,p3,labels=c("A","B"),ncol=1),p1,labels=c("","C"),ncol=2)
fig5
ggsave("Fig5.png",plot=fig5,height=16,width=20,units="cm",bg="white")


figS11<-ggplot(dat.tr.cor%>%filter(Trait=="gr_count")%>%
                 rename(c("PC1 climate"="PC1_climate","PC2 dimensions"="PC2_size",
                          "PC3 depth ratio"="PC3_depthRatio",
                          "PC4 latitude"="PC4_Latitude"))%>%
         pivot_longer(cols=starts_with("PC"),names_to="PC_axis",values_to="PC_val"),
       aes(x=PC_val,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("Precipitation seasonality")+ylab("Adjacent SMD\ngill rakers")+
  facet_wrap(vars(PC_axis))

figS11

ggsave("FigS11.png",plot=figS11,height=12,width=20,units="cm",bg="white")

figS12<-ggplot(dat.tr.cor%>%filter(Trait=="body_length")%>%
                 rename(c("PC1 climate"="PC1_climate","PC2 dimensions"="PC2_size",
                          "PC3 depth ratio"="PC3_depthRatio",
                          "PC4 latitude"="PC4_Latitude"))%>%
                 pivot_longer(cols=starts_with("PC"),names_to="PC_axis",values_to="PC_val"),
               aes(x=PC_val,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("Precipitation seasonality")+ylab("Adjacent SMD\nbody length")+
  facet_wrap(vars(PC_axis))

figS12

ggsave("FigS12.png",plot=figS12,height=12,width=20,units="cm",bg="white")

figS11<-ggplot(dat.tr.cor%>%filter(Trait=="gr_count")%>%
                 rename(c("PC1 climate"="PC1_climate","PC2 dimensions"="PC2_size",
                          "PC3 depth ratio"="PC3_depthRatio",
                          "PC4 latitude"="PC4_Latitude"))%>%
                 pivot_longer(cols=starts_with("PC"),names_to="PC_axis",values_to="PC_val"),
               aes(x=PC_val,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("Precipitation seasonality")+ylab("Adjacent SMD\ngill rakers")+
  facet_wrap(vars(PC_axis))

figS11

ggsave("FigS11.png",plot=figS11,height=12,width=20,units="cm",bg="white")

figS12<-ggplot(dat.tr.cor%>%filter(Trait=="body_length")%>%
                 rename(c("PC1 climate"="PC1_climate","PC2 dimensions"="PC2_size",
                          "PC3 depth ratio"="PC3_depthRatio",
                          "PC4 latitude"="PC4_Latitude"))%>%
                 pivot_longer(cols=starts_with("PC"),names_to="PC_axis",values_to="PC_val"),
               aes(x=PC_val,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("Precipitation seasonality")+ylab("Adjacent SMD\nbody length")+
  facet_wrap(vars(PC_axis))

figS12

ggsave("FigS12.png",plot=figS12,height=12,width=20,units="cm",bg="white")




p1<-ggplot(dat.tr.o.cor%>%filter(Trait=="gr_count"),
       aes(x=Prec_sea,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("Precipitation seasonality")+ylab("Outer SMD\ngill rakers")

p2<-ggplot(dat.tr.o.cor%>%filter(Trait=="gr_count"),
       aes(x=Temp_sea,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.cor))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("Temperature seasonality")+ylab("Outer SMD\ngill rakers")

p3<-ggplot(dat.tr.sal%>%filter(Trait=="body_length"),
       aes(x=Prec_sea,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.sal))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("Precipitation seasonality")+ylab("Adjacent SMD\nbody length")

p4<-ggplot(dat.tr.sal%>%filter(Trait=="body_length"),
       aes(x=logDepthMax,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.sal))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("Lake Depth")+ylab("Adjacent SMD\nbody length")

p5<-ggplot(dat.tr.o.sal%>%filter(Trait=="body_length"),
       aes(x=Prec_sea,y=yi,colour=Genus))+
  scale_colour_manual(values=c(col.sal))+
  geom_smooth(method="lm",alpha=0.2,size=2)+
  geom_point(size=3,shape=1,stroke=2,position=position_dodge(width=0.4))+
  theme(legend.position="none",text=element_text(size=12))+
  xlab("Precipitation seasonality")+ylab("Outer SMD\nbody length")
ggarrange(p1,p2,p3,p4,p5,nrow=3,ncol=2,labels=LETTERS)


morph_count.cv<-keepTraits%>%filter(Trait!="PC1_body")
morph_count.cv$slope.mean<-NA; morph_count.cv$slope.cI<-NA
morph_count.cv$slope.diff<-NA


for(i in 1:nrow(morph_count.cv)){
  lm.count<-lm(yi~Lacustrine_morphs,data=dat.tr.cv%>%
                 filter(Trait==morph_count.cv$Trait[i]&Genus==morph_count.cv$Genus[i]))
  morph_count.cv$slope.cI[i]<-paste(round(confint(lm.count)[2,1],3),
                                   round(confint(lm.count)[2,2],3),sep=", ")
  morph_count.cv$slope.mean[i]<-round(lm.count$coefficients[2],3)
  if(between(0,round(confint(lm.count)[2,1],3),
             round(confint(lm.count)[2,2],3))==F){
    morph_count.cv$slope.diff[i]<-"*" }else{
      morph_count.cv$slope.diff[i]<-"" }}

morph_count$Estimate<-"adjacent SMD";morph_count.o$Estimate<-"outer SMD";morph_count.cv$Estimate<-"lnCVR"
#write.csv(rbind(morph_count,morph_count.o,morph_count.cv),file="divergence by ecotype number.csv")


a<-ggplot(data=lm1.tr.o.weights,
       aes(x=pred_var,
           y=weight_sum,fill=genus))+
  facet_wrap(~genus+trait,ncol=2)+geom_bar(stat="identity")+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none",text=element_text(size=12))+
  ylim(0,1)+
  xlab("Abiotic factor")+ylab("Relative importance")+
  scale_x_discrete(labels=c("logDepthMax"="Depth","logArea"="Surface\narea",
                            "Temp_sea"="Temp\nseason","Prec_sea"="Precip\nseason"))+
  ggtitle("outer SMD")

b<-ggplot(data=lm1.tr.weights,
       aes(x=pred_var,
           y=weight_sum,fill=genus))+
  facet_wrap(~genus+trait,ncol=2)+geom_bar(stat="identity")+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none",text=element_text(size=12))+
  ylim(0,1)+
  xlab("Abiotic factor")+ylab("Relative importance")+
  scale_x_discrete(labels=c("logDepthMax"="Depth","logArea"="Surface\narea",
                            "Temp_sea"="Temp\nseason","Prec_sea"="Precip\nseason"))+
  ggtitle("adjacent SMD")

c<-ggplot(data=lm1.tr.cv.weights,
       aes(x=pred_var,
           y=weight_sum,fill=genus))+
  facet_wrap(~genus+trait,ncol=2)+geom_bar(stat="identity")+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none",text=element_text(size=12))+
  ylim(0,1)+
  xlab("Abiotic factor")+ylab("Relative importance")+
  scale_x_discrete(labels=c("logDepthMax"="Depth","logArea"="Surface\narea",
                            "Temp_sea"="Temp\nseason","Prec_sea"="Precip\nseason"))+
  ggtitle("lnCVR")


ggarrange(b,a,c,nrow=3)


#add "species" column to df
df<-df%>%left_join(m,by="Lake")%>%
  dplyr::select(Lake:Species)%>%distinct()
df$Genus<-substr(df$Species,1,3)

#remove NAs and divide dataframes by genus
df.cor<-df%>%filter(Genus=="Cor")%>%
  drop_na(all_of(c("PC1_climate","PC2_size","PC3_depthRatio","PC4_Latitude")))
df.sal<-df%>%filter(Genus=="Sal")%>%
  drop_na(all_of(c("PC1_climate","PC2_size","PC3_depthRatio","PC4_Latitude")))
  

#fit glm with quasi poisson because of underdispersion
lm2.cor<-glm(Lacustrine_morphs~PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
             family = quasipoisson(link="log"),
             data=df.cor,na.action="na.fail")
df.cor<-df.cor%>%
  mutate(ordinal_morphs=if_else(Lacustrine_morphs>4,4,Lacustrine_morphs))%>%
  mutate(ordinal_morphs=ordinal_morphs-1)


# brm.cor1<-brm(data = df.cor,
#     family = cumulative(probit),
#     ordinal_morphs ~ PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
#     prior(normal(0, 4), class = Intercept))

####Weights from weightable sum to 1
###Weight for each var is the sum of the weights of the models that include that var



df.sal<-df.sal%>%
  mutate(ordinal_morphs=if_else(Lacustrine_morphs>4,4,Lacustrine_morphs))%>%
  mutate(ordinal_morphs=ordinal_morphs-1)
unique(df.sal$ordinal_morphs)
# brm.sal1<-brm(data = df.sal,
#               family = cumulative(probit),
#               ordinal_morphs ~ PC1_climate+PC2_size+PC3_depthRatio+PC4_Latitude,
#               prior(normal(0, 4), class = Intercept))

#brm.cor<-fit_ordinal_models(df.cor)
#save(brm.cor, file = "CoregonusBrms.RData")
load("CoregonusBrms.RData")


#brm.sal<-fit_ordinal_models(df.sal)

#save(brm.sal, file = "SalvelinusBrms.RData")
load("SalvelinusBrms.RData")



#brm.list<-list(brm.cor,brm.cor2)
loo.cor<-list(loo(brm.cor[[1]]),loo(brm.cor[[2]]),loo(brm.cor[[3]]),
              loo(brm.cor[[4]]),loo(brm.cor[[5]]),loo(brm.cor[[6]]),
              loo(brm.cor[[7]]),loo(brm.cor[[8]]),loo(brm.cor[[9]]),
              loo(brm.cor[[10]]),loo(brm.cor[[11]]),loo(brm.cor[[12]]),
              loo(brm.cor[[13]]),loo(brm.cor[[14]]),loo(brm.cor[[15]]),
              loo(brm.cor[[16]]))

loo.sal<-list(loo(brm.sal[[1]]),loo(brm.sal[[2]]),loo(brm.sal[[3]]),
              loo(brm.sal[[4]]),loo(brm.sal[[5]]),loo(brm.sal[[6]]),
              loo(brm.sal[[7]]),loo(brm.sal[[8]]),loo(brm.sal[[9]]),
              loo(brm.sal[[10]]),loo(brm.sal[[11]]),loo(brm.sal[[12]]),
              loo(brm.sal[[13]]),loo(brm.sal[[14]]),loo(brm.sal[[15]]),
              loo(brm.sal[[16]]))

brm_factors<-read.table("pc_brm_inputs.txt")

calc_loo_weights<-function(looGenus,Genus.i,brm_factors){
  
  w1<-loo_model_weights(looGenus)
  
  
  colnames(brm_factors)<-c("model","factors")
  brm_factors$weights<-w1
  
  out.fl<-data.frame()
  for(pVar in tr_vars){ #cycle through model predictors
    #function uses grep, so check that no var names are subsets of other var names
    if(length(tr_vars[grepl(pVar,tr_vars)])>1){
      print("PROBLEM: variable name subset of another variable name")}
    df.temp<-data.frame(genus=Genus.i,trait="ecotypes", #save genus and trait name
                        pred_var=pVar, #save abiotic var name
                        #select all rows containing a trait name and sum those weights
                        AICc_weight=brm_factors%>%filter(grepl(pVar,factors))%>%
                          summarise(weight_sum=sum(weights)))
    out.fl<-rbind(out.fl,df.temp) #add results to dataframe
  }
  out.fl
}

lm2_weights<-rbind(calc_loo_weights(loo.cor,"Coregonus",brm_factors),
calc_loo_weights(loo.sal,"Salvelinus",brm_factors))
# write.csv(lm2_weights%>%mutate(weight_sum=round(weight_sum,3)),
#           "ecoCount.weights.240716.csv",row.names=F)


brm.sal.out<-brm(data = df.sal%>%rename("Sympatric_ecotypes"="ordinal_morphs"),
                 family = sratio("cloglog"),
                 Sympatric_ecotypes ~ PC1_climate+PC2_size)
summary(brm.sal.out)

col.cor.vals=c("#97C1C4",col.cor,"#0F2B2E")
col.sal.vals=c("#DAB3A7",col.sal,"#5A3427")

ce.sal.pc1 <- conditional_effects(brm.sal.out, "PC1_climate",categorical=T)
p1<-plot(ce.sal.pc1, plot = FALSE)[[1]] +
  scale_color_manual(values=col.sal.vals,
                     name="Sympatric\necotypes",
                   labels = c("2", "3", "4+")) +
  scale_fill_manual(values=col.sal.vals,
                  name="Sympatric\necotypes",
                  labels = c("2", "3", "4+"))+
  theme(text=element_text(size=12))+
  xlab("PC1 (climate)")

ce.sal.pc2 <- conditional_effects(brm.sal.out, "PC2_size",categorical=T)
p2<-plot(ce.sal.pc2, plot = FALSE)[[1]] +
  scale_color_manual(values=col.sal.vals,
                     name="Sympatric\necotypes",
                     labels = c("2", "3", "4+")) +
  scale_fill_manual(values=col.sal.vals,
                    name="Sympatric\necotypes",
                    labels = c("2", "3", "4+"))+
  theme(text=element_text(size=12))+
  xlab("PC2 (lake dimensions)")


brm.cor.out<-brm(data = df.cor,family = sratio("cloglog"),
                 ordinal_morphs ~ PC2_size+PC4_Latitude)
summary(brm.cor.out)

ce.cor.pc2 <- conditional_effects(brm.cor.out, "PC2_size",categorical=T)
p3<-plot(ce.cor.pc2, plot = FALSE)[[1]] +
  scale_color_manual(values=col.cor.vals,
                     name="Sympatric\necotypes",
                     labels = c("2", "3", "4+")) +
  scale_fill_manual(values=col.cor.vals,
                    name="Sympatric\necotypes",
                    labels = c("2", "3", "4+"))+
  theme(text=element_text(size=12))+
  xlab("PC2 (lake dimensions)")

ce.cor.pc2 <- conditional_effects(brm.cor.out, "PC4_Latitude",categorical=T)
p4<-plot(ce.cor.pc2, plot = FALSE)[[1]] +
  scale_color_manual(values=col.cor.vals,
                     name="Sympatric\necotypes",
                     labels = c("2", "3", "4+")) +
  scale_fill_manual(values=col.cor.vals,
                    name="Sympatric\necotypes",
                    labels = c("2", "3", "4+"))+
  theme(text=element_text(size=12))+
  xlab("PC4 (latitude)")

ggarrange(ggarrange(p1,p2,common.legend=T,legend="right",labels=c("A","B")),
          ggarrange(p4,p3,common.legend=T,legend="right",labels=c("C","D")),
          nrow=2,ncol=1)

p5<-ggplot(data=lm2_weights,aes(x=pred_var,
                            y=weight_sum,fill=genus))+
  facet_wrap(~genus,ncol=1)+geom_bar(stat="identity")+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none",text=element_text(size=12))+
  ylim(0,1)+
  xlab("Abiotic factor")+ylab("Relative importance")+
  scale_x_discrete(labels=c("PC1_climate"="PC1\nclimate","PC2_size"="PC2\ndimensions",
                            "PC3_depthRatio"="PC3\ndepth ratio",
                            "PC4_Latitude"="PC4\nlatitude"))


ggarrange(ggarrange(ggarrange(p3,p4,common.legend=T,legend="right",labels=c("A","")),
  ggarrange(p2,p1,common.legend=T,legend="right",labels=c("B","")),
          nrow=2,ncol=1),
  p5,nrow=1,ncol=2,widths=c(1.5,1),labels=c("","C"))




#plotting results
p2<-ggplot(data=lm2_weights,aes(x=reorder(pred_var,-weight_sum),
                            y=weight_sum,fill=genus))+
  facet_wrap(~genus,ncol=1)+geom_bar(stat="identity")+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(legend.position="none",text=element_text(size=12))+
  ylim(0,1)+
  xlab("Abiotic factor")+ylab("Relative importance")+
  scale_x_discrete(labels=c("logDepthMax"="Depth","logArea"="Surface\narea",
                            "Temp_mean"="Temp","Temp_sea"="Temp\nseason",
                            "Prec_mean"="Precip","Prec_sea"="Precip\nseason",
                            "Lat"="Latitude"))

p1<-ggplot(data=df,aes(x=logArea,y=logDepthMax,size=Lacustrine_morphs,colour=Genus))+
  geom_point()+
  scale_colour_manual(values=c(col.cor,col.sal))+
  theme(text=element_text(size=12),legend.position.inside=c(0.85,0.2))+
  xlab("Surface Area (log scale)")+ylab("Depth (log scale)")+
  labs(size="Ecotypes")+guides(colour="none")

p1a<-ggplot(data=df,aes(x=PC2_size,y=Lacustrine_morphs,size=Lacustrine_morphs,fill=Genus,colour=Genus))+
  #geom_smooth(method="lm")+
  geom_point()+
  scale_colour_manual(values=c(col.cor,col.sal))+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(text=element_text(size=12),legend.position="none")+
  xlab("PC2 (lake size)")+ylab("Sympatric ecotypes")+
  labs(size="Ecotypes")+guides(colour="none")
p1b<-ggplot(data=df,aes(x=PC1_climate,y=Lacustrine_morphs,size=Lacustrine_morphs,fill=Genus,colour=Genus))+
  #geom_smooth(method="lm")+
  geom_point()+
  scale_colour_manual(values=c(col.cor,col.sal))+
  scale_fill_manual(values=c(col.cor,col.sal))+
  theme(text=element_text(size=12),legend.position="none")+
  xlab("PC1 (climate)")+ylab("Sympatric ecotypes")+
  labs(size="Ecotypes")+guides(colour="none")

ggarrange(ggarrange(p1b,p1a,labels=LETTERS,nrow=2,ncol=1),p2,labels=c("","C"))


###calculate mean traits

mean_traits<-as.data.frame(t(rep("row1",7))) #adjacent dataframe
colnames(mean_traits)<-cbind("mean","CI.l","CI.u","trait","B_nEco","CI.l_nEco","CI.u_nEco")
mean_traits.o<-as.data.frame(t(rep("row1",7))) #outer dataframe
colnames(mean_traits.o)<-cbind("mean","CI.l","CI.u","trait","B_nEco","CI.l_nEco","CI.u_nEco")
mean_traits.cv<-as.data.frame(t(rep("row1",7))) #lnCVR dataframe
colnames(mean_traits.cv)<-cbind("mean","CI.l","CI.u","trait","B_nEco","CI.l_nEco","CI.u_nEco")

mean_traits<-add_to_meantraits(dat.gr,"Gill\nRakers",mean_traits,F) #gill rakers adjacent
mean_traits.o<-add_to_meantraits(dat.gr.o,"Gill\nRakers",mean_traits.o,F) #gill rakers outer
mean_traits.cv<-add_to_meantraits(dat.gr.cv,"Gill\nRakers",mean_traits.cv,F) #gill rakers outer

mean_traits<-add_to_meantraits(dat.fl,"Body\nLength",mean_traits,F) #body length adjacent
mean_traits.o<-add_to_meantraits(dat.fl.o,"Body\nLength",mean_traits.o,F) #body length outer
mean_traits.cv<-add_to_meantraits(dat.fl.cv,"Body\nLength",mean_traits.cv,F) #body length outer

mean_traits<-add_to_meantraits(dat.ag,"Age",mean_traits,F) #age adjacent
mean_traits.o<-add_to_meantraits(dat.ag.o,"Age",mean_traits.o,F) #age outer
mean_traits.cv<-add_to_meantraits(dat.ag.cv,"Age",mean_traits.cv,F) #age outer

mean_traits<-add_to_meantraits(dat.al,"Asymptotic\nLength",mean_traits,F) #Asymptotic adjacent
mean_traits.o<-add_to_meantraits(dat.al.o,"Asymptotic\nLength",mean_traits.o,F) #Asymptotic length outer
mean_traits.cv<-add_to_meantraits(dat.al.cv,"Asymptotic\nLength",mean_traits.cv,F) #Asymptotic length outer

mean_traits<-add_to_meantraits(dat.N15,"Isotope\nN15",mean_traits,F) #N15 adjacent
mean_traits.o<-add_to_meantraits(dat.N15.o,"Isotope\nN15",mean_traits.o,F) #N15 outer
mean_traits.cv<-add_to_meantraits(dat.N15.cv,"Isotope\nN15",mean_traits.cv,F) #N15 lnCVR

mean_traits<-add_to_meantraits(dat.C13,"Isotope\nC13",mean_traits,F) #C13 adjacent
mean_traits.o<-add_to_meantraits(dat.C13.o,"Isotope\nC13",mean_traits.o,F) #C13 outer
mean_traits.cv<-add_to_meantraits(dat.C13.cv,"Isotope\nC13",mean_traits.cv,F) #C13 lnCVR

mean_traits<-add_to_meantraits(dat.troph,"Head\nLengths",mean_traits,T) #trophic adjacent
mean_traits.o<-add_to_meantraits(dat.troph.o,"Head\nLengths",mean_traits.o,T) #trophic outer
mean_traits.cv<-add_to_meantraits(dat.troph.cv,"Head\nLengths",mean_traits.cv,T) #trophic lnCVR

mean_traits<-add_to_meantraits(dat.fin,"Fin\nLengths",mean_traits,T) #fins adjacent
mean_traits.o<-add_to_meantraits(dat.fin.o,"Fin\nLengths",mean_traits.o,T) #fins outer
mean_traits.cv<-add_to_meantraits(dat.fin.cv,"Fin\nLengths",mean_traits.cv,T) #fins lnCVR

mean_traits<-add_to_meantraits(dat.dep,"Capture\nDepth",mean_traits,F) #depth adjacent
mean_traits.o<-add_to_meantraits(dat.dep.o,"Capture\nDepth",mean_traits.o,F) #depth outer
mean_traits.cv<-add_to_meantraits(dat.dep.cv,"Capture\nDepth",mean_traits.cv,F) #depth lnCVR

mean_traits<-add_to_meantraits(dat.pcb,"PCA-body\nshape",mean_traits,F) #PCA body adjacent
mean_traits.o<-add_to_meantraits(dat.pcb.o,"PCA-body\nshape",mean_traits.o,F) #PCA body outer

mean_traits<-add_to_meantraits(dat.pch,"PCA-head\nshape",mean_traits,F) #PCA head adjacent
mean_traits.o<-add_to_meantraits(dat.pch.o,"PCA-head\nshape",mean_traits.o,F) #PCA head outer

#set necessary columns as numeric
mean_traits[,1:3]<-mutate_all(mean_traits[,1:3],as.numeric)
mean_traits$Genus<-t(data.frame(strsplit(mean_traits$trait," ")))[,2]
mean_traits$trait<-t(data.frame(strsplit(mean_traits$trait," ")))[,1]

#graph mean trait effect sizes

p1<-ggplot(data=mean_traits,aes(x=trait,y=mean,ymin=CI.l,ymax=CI.u,colour=Genus,shape=Genus))+
  geom_point(size=6,position=position_dodge(0.4))+
  geom_errorbar(position=position_dodge(0.4),width=0.3,size=1.5)+
  ylab("Adjacent SMD")+xlab("Trait")+
  scale_colour_manual(values=c(col.cor,col.sal))+
  scale_shape_manual(values=c(16,17))+
  theme(text=element_text(size=12))+
  theme(legend.position=c(.9,.8))
p1

mean_traits%>%group_by(trait)%>%dplyr::select(trait,mean,Genus)%>%
  pivot_wider(names_from = Genus,values_from = mean)%>%
  mutate(Cor.diff=Coregonus/Salvelinus,Sal.diff=Salvelinus/Coregonus)


#set necessary columns as numeric
mean_traits.o[,1:3]<-mutate_all(mean_traits.o[,1:3],as.numeric)
mean_traits.o$Genus<-t(data.frame(strsplit(mean_traits.o$trait," ")))[,2]
mean_traits.o$trait<-t(data.frame(strsplit(mean_traits.o$trait," ")))[,1]

mean_traits.o%>%group_by(trait)%>%dplyr::select(trait,mean,Genus)%>%
  pivot_wider(names_from = Genus,values_from = mean)%>%
  mutate(Cor.diff=Coregonus/Salvelinus,Sal.diff=Salvelinus/Coregonus)


p2<-ggplot(data=mean_traits.o,aes(x=trait,y=mean,ymin=CI.l,ymax=CI.u,colour=Genus,shape=Genus))+
  geom_point(size=6,position=position_dodge(0.4))+
  geom_errorbar(position=position_dodge(0.4),width=0.3,size=1.5)+
  ylab("outer SMD")+xlab("Trait")+
  scale_colour_manual(values=c(col.cor,col.sal))+
  scale_shape_manual(values=c(16,17))+
  theme(text=element_text(size=12))+
  theme(legend.position="none")
p2

mean_traits.cv[,1:3]<-mutate_all(mean_traits.cv[,1:3],as.numeric)
mean_traits.cv$Genus<-t(data.frame(strsplit(mean_traits.cv$trait," ")))[,2]
mean_traits.cv$trait<-t(data.frame(strsplit(mean_traits.cv$trait," ")))[,1]

p3<-ggplot(data=mean_traits.cv,aes(x=trait,y=mean,ymin=CI.l,ymax=CI.u,colour=Genus,shape=Genus))+
  geom_point(size=6,position=position_dodge(0.4))+
  geom_errorbar(position=position_dodge(0.4),width=0.3,size=1.5)+
  ylab("lnCVR")+xlab("Trait")+
  scale_colour_manual(values=c(col.cor,col.sal))+
  theme(text=element_text(size=12))+
  scale_shape_manual(values=c(16,17))

p3
ggarrange(p2,p3,labels=c("A","B"),nrow=2)

lakes.trOnly<-dat.tr%>%pull(Lake)%>%unique()

ggplot(data=dat.tr,aes(y=yi,x=Lacustrine_morphs,colour=Genus,shape=Genus))+
  facet_wrap(~Trait,nrow=2)+
  geom_point(size=3,position=position_dodge(0.4))+
  geom_smooth(method="lm",alpha=0)+
  scale_colour_manual(values=c(col.cor,col.sal))+
  theme(text=element_text(size=12))+
  scale_shape_manual(values=c(16,17))+
  xlab("Number of ecotypes")+ylab("Adjacent SMD")

ggplot(data=dat.tr.o,aes(y=yi,x=Lacustrine_morphs,colour=Genus,shape=Genus))+
  facet_wrap(~Trait,nrow=2)+
  geom_point(size=3,position=position_dodge(0.4))+
  geom_smooth(method="lm",alpha=0)+
  scale_colour_manual(values=c(col.cor,col.sal))+
  theme(text=element_text(size=12))+
  scale_shape_manual(values=c(16,17))+
  xlab("Number of ecotypes")+ylab("Outer SMD")

ggplot(data=dat.tr.cv,aes(y=yi,x=Lacustrine_morphs,colour=Genus,shape=Genus))+
  facet_wrap(~Trait,nrow=2)+
  geom_point(size=3,position=position_dodge(0.4))+
  geom_smooth(method="lm",alpha=0)+
  scale_colour_manual(values=c(col.cor,col.sal))+
  theme(text=element_text(size=12))+
  scale_shape_manual(values=c(16,17))+
  xlab("Number of ecotypes")+ylab("lnCVR")
  


#count up the lakes in the study
lakes.n<-df1$Lake%>%unique()

lakes.all<-c(lakes.n,lakes.tr,lakes.trOnly)%>%unique()

