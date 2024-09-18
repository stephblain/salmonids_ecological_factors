#Functions

#calculates mean and variance for pairs
##uses escalc function which is equivalent to eqns 1 and 2 from Lajeunesse 2011
##Input: dataframe for one study/lake combination, with one row per morph,
###and columns named "Lake", "mean", "n", "sd", and "group"
calc.y<-function(x){co<-cbind(1:(nrow(x)-1),2:nrow(x)) #get matrix to order everything else
x<-x[order(-x$mean),] #order by trait mean to get correct comparisons
x1<-data.frame(m1i=x$mean[co[,1]],m2i=x$mean[co[,2]], #make matrix with means, sd's and sample sizes
               n1i=x$n[co[,1]],n2i=x$n[co[,2]],sd1i=x$sd[co[,1]],sd2i=x$sd[co[,2]],
               group=x$group[co[,1]])
#Remove morph if sample size is less than 5
x1<-x1[x1$n1i>=5,];x1<-x1[x1$n2i>=5,]
#print error message if no comparison to morph with sample size greater than 5
if(length(x1$m1i)==0){print(paste("sample size too small for",x$group[1]))}else{
  x1<-cbind(x1,as.data.frame(escalc(measure="SMD",m1i=x1$m1i,m2i=x1$m2i, #use escalc for SMD and var calculations
                                    n1i=x1$n1i,n2i=x1$n2i, ## returns the same result as manual eqns below
                                    sd1i=x1$sd1i,sd2i=x1$sd2i)))
  attr(x1$yi,"measure")<-NULL;attr(x1$yi,"ni")<-NULL #remove attributes
  x1}}
# #Manual calculation with bias correction based on d in Hedges and Olkin 1985
# sd_pooled<-sqrt(((x1$n1i-1)*(x1$sd1i^2)+(x1$n2i-1)*(x1$sd2i^2))/(x1$n1i+x1$n2i-2))
# N<-x1$n1i+x1$n2i
# g<-((x1$m1i-x1$m2i)/sd_pooled)
# J<-1-3/(4*(N-2)-1)
# yi<-J*g
# #Variance  for d from hedges 1985 and gleser 2009
# N/(x1$n1i*x1$n2i)+(yi^2)/(2*N) == 1/x1$n1i + 1/x1$n2i+yi^2/(2*N)




#Calculates variance/covariance matrix using eqns 19.18 (variance) and 19.19 (covariance) from Gleser and Olkin 2009
##the difference is that morph pairs with no overlapping individuals are given
##a covariance of zero
#"control" covariance is used between morph pairs that share a morph
##Input: output from calc.y
calc.v<-function(x){n<-nrow(x) # number of effect sizes
v<-matrix(0,nrow=n, ncol=n) # zero matrix
N<-sum(x$n1i,x$n2i[n]) #total sample size (all n1's plus the last n2 = all ecotypes)
#calculate variances with total N (if n1 + n2 used in place of N, matches escalc output)
x$vi<-1/x$n1i + 1/x$n2i + x$yi^2/(2*N) #from eqn 19.18 Gleser and Olkin 2009
for(i in 1:(n-1)){ # calculate and insert covariances
  v1<-1/x$n2i[i]+x$yi[i]*x$yi[(i+1)]/(2*N) #from eqn 19.19 Gleser and Olkin 2009
  v[i,(i+1)]<-v1;v[(i+1),i]<-v1} #insert covariances
diag(v)<-x$vi;v} # insert variances

#outputs an aggregated mean and variance for multi-morph lakes
#uses eqn 3 from Lajeunesse 2011
#input: same as calc.y
calc.multi<-function(x){
  if(is.numeric(calc.y(x)$yi)){
    V<-calc.v(calc.y(x))
    X<-t(t(rep(1,nrow(V))))
    E<-t(t(calc.y(x)$yi))
    if(is.positive.definite(V)){
      SMD<-solve(t(X)%*%solve(V)%*%X)*(t(X)%*%solve(V)%*%E) #eqn 19.21 from Gleser and Olkin 2009
      #(and eqn 3 from Lajeunesse 2008)
      #calculation has the same result as rma.mv
      vSMD<-solve(t(X)%*%solve(V)%*%X)}else{ #also eqn 19.21 (and paragraph below eqn 3)
        SMD<-NA
        vSMD<-NA}
    out<-cbind(x$group[1],SMD,vSMD); colnames(out)<-c("group","yi","vi");out}else{print(calc.y(x))}}

#outputs a dataframe with an effect size and variance for each lake, for a trait
#input: column numbers from morphs dataframe in order (1) mean, (2) n, (3) e, (4) e type
#       traitname - character string
calc.dat.tr<-function(tr_cols,traitname){


  tr<-morphs[c(1:6,tr_cols)]
  
  colnames(tr)[7:10]<-c("tr.mean","tr.n","tr.e","tr.etype")
  tr$tr.mean<-as.numeric(tr$tr.mean)
  tr$tr.n<-as.numeric(tr$tr.n)
  tr$tr.e<-as.numeric(tr$tr.e)
  tr<-tr[!is.na(tr$tr.n),] # remove if no sample size
  tr<-tr[!is.na(tr$tr.e),] # remove if no error
  #convert Variance, SE and 95CI to SD
  tr$tr.e[tr$tr.etype=="Variance"]<-sqrt(tr$tr.e[tr$tr.etype=="Variance"])
  tr$tr.e[tr$tr.etype=="SE"]<-se_to_sd(tr$tr.e[tr$tr.etype=="SE"],
                                       tr$tr.n[tr$tr.etype=="SE"])
  tr$tr.e[tr$tr.etype=="95CI"]<-CI95_to_sd(tr$tr.e[tr$tr.etype=="95CI"],
                                           tr$tr.n[tr$tr.etype=="95CI"])
  tr$group<-paste(tr$Paper_ID,tr$Lake,sep="_") #make group factor, unique to each lake/study combo
  groups<-levels(as.factor(tr$group))
  dat.tr<-data.frame("row1","row1","row1");colnames(dat.tr)<-c("group","yi","vi")
  for(i in groups){z<-tr[tr$group==i,c(2,7:9,11)] #select needed columns
  names(z)<-c("Lake","mean","n","sd","group") #rename columns
  #check that there are two ecotypes with sufficient data in the lake
  if(nrow(z)<2){print(paste("sample size too small for",i))}else{
    calc.y.z<-calc.y(z) #calculated SMD + variance for each adjacent pair
    if(is.character(calc.y.z)){print(calc.y.z)}else{
      if(nrow(calc.y.z)==1){ # if 2-ecotype, add to final dataframe
        out<-as.matrix(calc.y.z[,7:9])
        dat.tr<-rbind(dat.tr,out)}
      if(nrow(calc.y.z)>1){ #if 3+ ecotype, run through calc.multi to get estimate for lake
        out<-calc.multi(z)
        dat.tr<-rbind(dat.tr,out)}}}}
  rm(out,z,i)
  dat.tr<-dat.tr[-1,] #remove "row1" row
  dat.tr$group<-as.character(dat.tr$group)
  dat.tr$yi<-as.numeric(as.character(dat.tr$yi))
  dat.tr$vi<-as.numeric(as.character(dat.tr$vi))
  dat.tr$Lake<-t(data.frame(strsplit(dat.tr$group,"_")))[,2]
  
  #fixed effects meta-analysis to get one estimate per lake
  tab1<-data.frame(table(dat.tr$Lake)) #counts of estimates
  tab1<-tab1[tab1$Freq>1,] #get lakes with repeated estimates
  if(nrow(tab1)>0){ #if there are multiple estimates for a lake
    for(x in 1:nrow(tab1)){ #run through repeated estimate lakes
      dat.sub1<-dat.tr[dat.tr$Lake==tab1$Var1[x],]
      res.sub1<-rma(dat.sub1$yi,dat.sub1$vi,method="FE") #fixed effects meta-analysis
      #results formatting to add to dat.tr
      Lake<-dat.sub1$Lake[1]
      dattr.sub1<-data.frame(cbind(paste("01E",Lake,sep="_"),res.sub1$beta,
                                   res.sub1$vb, #sqrt(vi)==res.sub1$se
                                   Lake))
      colnames(dattr.sub1)<-colnames(dat.tr); rownames(dattr.sub1)<-NULL
      dattr.sub1$yi<-as.numeric(dattr.sub1$yi); dattr.sub1$vi<-as.numeric(dattr.sub1$vi)
      dat.tr<-dat.tr[!dat.tr$Lake==tab1$Var1[x],] #remove lake from dat.tr
      dat.tr<-rbind(dat.tr,dattr.sub1)}} #add new estimate to dat.tr
  
  #add lake data to dat.tr dataframe
  dat.tr<-merge(dat.tr,df,by="Lake",all.x=T)
  dat.tr<-dat.tr[,c("Lake","group","yi","vi","Country","River_system",
                    "Lat","Lon","logArea","logDepthMax","Lacustrine_morphs",
                    "Temp_mean","Temp_sea","Prec_mean","Prec_sea",
                    "Lat.2","logArea.2","logDepthMax.2",
                    "Temp.2","Prec.2","Continent",
                    "PC1_climate","PC2_size","PC3_depthRatio","PC4_Latitude")]
  dat.tr<-merge(dat.tr,distinct(tr[,c(2,3)]),by="Lake",all.y=T)
  if(nrow(dat.tr)==length(unique(tr$Lake))){print("merging probably worked okay")}else{
    print("something went wrong, check merging")}
  dat.tr<-dat.tr[!is.na(dat.tr$vi),]
  dat.tr$sdi<-sqrt(dat.tr$vi)
  dat.tr$Genus<-t(data.frame(strsplit(dat.tr$Species,"_")))[,1]
  dat.tr$Trait<-traitname; dat.tr}


#function to calculate the effect size using only the most extreme morphs in a lake
#i.e. estimate effect size of the breadth of divergence
calc.dat.tr_outer<-function(tr_cols,traitname){
  
  
  tr<-morphs[c(1:4,6:7,tr_cols)]
  colnames(tr)[7:10]<-c("tr.mean","tr.n","tr.e","tr.etype")
  tr$tr.mean<-as.numeric(tr$tr.mean)
  tr$tr.n<-as.numeric(tr$tr.n)
  tr$tr.e<-as.numeric(tr$tr.e)
  tr<-tr[!is.na(tr$tr.n),] # remove if no sample size
  tr<-tr[!is.na(tr$tr.e),] # remove if no error
  #convert Variance, SE and 95CI to SD
  tr$tr.e[tr$tr.etype=="Variance"]<-sqrt(tr$tr.e[tr$tr.etype=="Variance"])
  tr$tr.e[tr$tr.etype=="SE"]<-se_to_sd(tr$tr.e[tr$tr.etype=="SE"],
                                       tr$tr.n[tr$tr.etype=="SE"])
  tr$tr.e[tr$tr.etype=="95CI"]<-CI95_to_sd(tr$tr.e[tr$tr.etype=="95CI"],
                                           tr$tr.n[tr$tr.etype=="95CI"])
  tr$group<-paste(tr$Paper_ID,tr$Lake,sep="_") #make group factor, unique to each lake/study combo
  groups<-levels(as.factor(tr$group))
  dat.tr<-data.frame("row1","row1","row1");colnames(dat.tr)<-c("group","yi","vi")
  
  for(i in groups){
    z<-tr[tr$group==i,c(2,7:9,11)] #select needed columns
    names(z)<-c("Lake","mean","n","sd","group") #rename columns
    if(nrow(z)>1){
      calc.y.z<-calc.y(z)
      if(is.character(calc.y.z)){print(calc.y.z)}else{
        if(nrow(calc.y.z)==1){
          out<-as.matrix(calc.y.z[,7:9])
          dat.tr<-rbind(dat.tr,out)}
        if(nrow(calc.y.z)>1){
          z<-z[order(z$mean),];z<-z[z$n>=5,]
          z<-z[c(1,nrow(z)),]
          out<-as.matrix(calc.y(z)[,7:9])
          dat.tr<-rbind(dat.tr,out)}}}}
  rm(out,z,i)
  dat.tr<-dat.tr[-1,]
  dat.tr$group<-as.character(dat.tr$group)
  dat.tr$yi<-as.numeric(as.character(dat.tr$yi))
  dat.tr$vi<-as.numeric(as.character(dat.tr$vi))
  dat.tr$Lake<-t(data.frame(strsplit(dat.tr$group,"_")))[,2]
  
  #fixed effects meta-analysis to get one estimate per lake
  tab1<-data.frame(table(dat.tr$Lake)) #counts of estimates
  tab1<-tab1[tab1$Freq>1,] #get lakes with repeated estimates
  if(nrow(tab1)>0){ #if there are multiple estimates for a lake
    for(x in 1:nrow(tab1)){ #run through repeated estimate lakes
      dat.sub1<-dat.tr[dat.tr$Lake==tab1$Var1[x],]
      res.sub1<-rma(dat.sub1$yi,dat.sub1$vi,method="FE") #fixed effects meta-analysis
      #results formatting to add to dat.tr
      Lake<-dat.sub1$Lake[1]
      dattr.sub1<-data.frame(cbind(paste("01E",Lake,sep="_"),res.sub1$beta,
                                   res.sub1$vb, #sqrt(vi)==res.sub1$se
                                   Lake))
      colnames(dattr.sub1)<-colnames(dat.tr); rownames(dattr.sub1)<-NULL
      dattr.sub1$yi<-as.numeric(dattr.sub1$yi); dattr.sub1$vi<-as.numeric(dattr.sub1$vi)
      dat.tr<-dat.tr[!dat.tr$Lake==tab1$Var1[x],] #remove lake from dat.tr
      dat.tr<-rbind(dat.tr,dattr.sub1)}} #add new estimate to dat.tr
  
  #add lake data to dat.gr dataframe
  dat.tr<-merge(dat.tr,df,by="Lake",all.x=T)
  dat.tr<-dat.tr[,c("Lake","group","yi","vi","Country","River_system",
                    "Lat","Lon","logArea","logDepthMax","Lacustrine_morphs",
                    "Temp_mean","Temp_sea","Prec_mean","Prec_sea",
                    "Lat.2","logArea.2","logDepthMax.2",
                    "Temp.2","Prec.2","Continent",
                    "PC1_climate","PC2_size","PC3_depthRatio","PC4_Latitude")]
  dat.tr<-merge(dat.tr,distinct(tr[,c(2,3)]),by="Lake",all.y=T)
  if(nrow(dat.tr)==length(unique(tr$Lake))){print("merging probably worked okay")}else{
    print("something went wrong, check merging")}
  dat.tr<-dat.tr[!is.na(dat.tr$yi),]#remove lakes without yi calc
  dat.tr$sdi<-sqrt(dat.tr$vi)
  dat.tr$Genus<-t(data.frame(strsplit(dat.tr$Species,"_")))[,1]
  dat.tr$Trait<-traitname; dat.tr}



#function to make graphs relating trait values to lake summary data
#input: output from calc.multi and trait name
graphs_lakeVars<-function(dat.tr,trait_name){
  species_colour<-data.frame(cbind(levels(as.factor(morphs$Species)),
                                   c(brewer.pal(8,"BuGn")[c(5,6,7,8,4)],brewer.pal(8,"OrRd")[c(8,7,6,5)])))
  tr_colour<-as.character(species_colour[species_colour$X1%in%levels(as.factor(dat.tr$Species)),2])
  p1<-ggplot(data=dat.tr,aes(x=Lacustrine_morphs,y=yi,colour=Species))+
    geom_pointrange(aes(ymin=(yi-sdi),ymax=(yi+sdi)),position=position_jitter(0.1),pch=19,size=0.4)+
    theme_classic()+
    scale_color_manual(values=tr_colour)+
    ylab(trait_name)+xlab("Number of Morphs")
  p2<-ggplot(data=dat.tr,aes(x=Lat,y=yi,colour=Species))+
    geom_pointrange(aes(ymin=(yi-sdi),ymax=(yi+sdi)),position=position_jitter(0.1),pch=19,size=0.4)+
    theme_classic()+
    scale_color_manual(values=tr_colour)+
    ylab(trait_name)+xlab("Latitude")
  p3<-ggplot(data=dat.tr,aes(x=logArea,y=yi,colour=Species))+
    geom_pointrange(aes(ymin=(yi-sdi),ymax=(yi+sdi)),position=position_jitter(0.1),pch=19,size=0.4)+
    theme_classic()+
    scale_color_manual(values=tr_colour)+
    ylab(trait_name)+xlab("ln Area")
  p4<-ggplot(data=dat.tr,aes(x=logDepthMax,y=yi,colour=Species))+
    geom_pointrange(aes(ymin=(yi-sdi),ymax=(yi+sdi)),position=position_jitter(0.1),pch=19,size=0.4)+
    theme_classic()+
    scale_color_manual(values=tr_colour)+
    ylab(trait_name)+xlab("ln Depth")
  grid.arrange(p1,p2,p3,p4,nrow=2)}


#dat.tr is the dataframe outputted from calc.dat.tr
#traitname is the name of the trait being added
#phylo is a Newick format phylogeny
#meantraits is a preexisting 4 column dataframe with column names/order "mean","CI.l","CI.u","trait"
#meantraits should have some sort of first row, it will be removed if mean is "row1"
#this function uses rm.mv to do a multivariate meta-analysis of the data for a trait
##with Lake and phylo as random effects if multi=F
##with group, Lake, and phylo as random effects if multi=T
add_to_meantraits<-function(dat.tr,traitname,meantraits,multi){

  ##Temporarily un-comment for editing code
  # dat.tr<-dat.fl
  # traitname<-"Body\nLength"
  # meantraits<-mean_traits
  # multi<-F

  #make column numeric if necessary
  dat.tr$Lacustrine_morphs<-as.numeric(as.character(dat.tr$Lacustrine_morphs))
  
  dat.tr<-dat.tr[!is.na(dat.tr$Lake),]
  #phylo.corr<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,dat.tr$Species))))
  dat.tr.C<-subset(dat.tr,dat.tr$Genus=="Coregonus")
  # if(nrow(dat.tr.C)>5){
  #   phylo.corr.C<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,dat.tr.C$Species))))}
  # 
  dat.tr.S<-subset(dat.tr,dat.tr$Genus=="Salvelinus")
  # if(nrow(dat.tr.S>5)){
  #   phylo.corr.S<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,dat.tr.S$Species))))}
  # 
  #start for multiple traits
  if(multi==T){
    if(nrow(dat.tr.C)>5){rma_tr.C<-rma.mv(yi,vi,random=list(~1|Lake,~1|Trait),
                                          data=dat.tr.C,method="REML")}
    if(nrow(dat.tr.S)>5){rma_tr.S<-rma.mv(yi,vi,random=list(~1|Lake,~1|Trait),
                                        data=dat.tr.S,method="REML")}
    
    if(nrow(dat.tr.C)>5){rma_nEco.C<-rma.mv(yi~Lacustrine_morphs,vi,
                                            random=list(~1|Lake,~1|Trait),
                                          data=dat.tr.C,method="REML")}
    if(nrow(dat.tr.S)>5){rma_nEco.S<-rma.mv(yi~Lacustrine_morphs,vi,
                                            random=list(~1|Lake,~1|Trait),
                                          data=dat.tr.S,method="REML")}}
  #start for one trait only
  if(multi==F){
    if(nrow(dat.tr.C)>5){rma_tr.C<-rma.uni(yi,vi,
                                          data=dat.tr.C,method="REML")}
    if(nrow(dat.tr.S)>5){rma_tr.S<-rma.uni(yi,vi,
                                          data=dat.tr.S,method="REML")}
    
    if(nrow(dat.tr.C)>5){rma_nEco.C<-rma.uni(yi~Lacustrine_morphs,vi,
                                           data=dat.tr.C,method="REML")}
    if(nrow(dat.tr.S)>5){rma_nEco.S<-rma.uni(yi~Lacustrine_morphs,vi,
                                           data=dat.tr.S,method="REML")}}
  
  
  
  meantraits<-mutate_all(meantraits,as.character)
  
  if(nrow(dat.tr.C)>5){
    meantraits<-rbind(meantraits,c(rma_tr.C$beta,rma_tr.C$ci.lb,rma_tr.C$ci.ub,
                                   paste(traitname,"Coregonus"),
                                   rma_nEco.C$beta[2],rma_nEco.C$ci.lb[2],rma_nEco.C$ci.ub[2]))}
  
  if(nrow(dat.tr.S)>5){
    meantraits<-rbind(meantraits,c(rma_tr.S$beta,rma_tr.S$ci.lb,rma_tr.S$ci.ub,
                                   paste(traitname,"Salvelinus"),
                                   rma_nEco.S$beta[2],rma_nEco.S$ci.lb[2],rma_nEco.S$ci.ub[2]))}
  
  if(meantraits$mean[1]=="row1"){meantraits<-meantraits[-1,]}
  meantraits}




#Data! Analysis!
#Input:
#dat.tr is the dataframe outputted from calc.dat.tr
#EStype is the type of effect size, generally "adjacent" or "outer"
#traitname is the name of the trait being added
#phylo is a Newick format phylogeny with names that correspond to Species column of dat.tr
#print.aic is T or F and determines whether to output AIC or anova results
#Output: AIC comparison table or Anova results of a random effects multivariate met analysis
calc.res.tr<-function(dat.tr,EStype,traitname,phylo,print.aic){

  dat.tr<-dat.tr[!is.na(dat.tr$yi),]
  res.tr<-data.frame("row1","row1","row1","row1","row1","row1","row1")
  res.aic<-data.frame("row1","row1","row1","row1","row1","row1")
  #drop unnecessary tips, compute branch lengths, make correlation matrix
  phylo.corr<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,dat.tr$Species))))
  colnames(res.aic)<-c("trait","feature","ES_type","model","AIC","AIC.diff")
  colnames(res.tr)<-c("trait","feature","pred","ES_type","estimate","CI.lower","CI.upper")
  for(x in c("Lat","logDepthMax","logArea","Lacustrine_morphs")){
    
    dat.tr.x<-dat.tr[,x]
    dat.tr.x2<-dat.tr.x^2
    #ML for model selection
    temp.1<-rma.mv(yi,vi,mods=~dat.tr.x+dat.tr.x2,random=list(~1|Species),
                   R=list(Species=phylo.corr),
                   data=dat.tr,method="ML")
    
    temp.2<-rma.mv(yi,vi,mods=~dat.tr.x,random=list(~1|Species),
                   R=list(Species=phylo.corr),
                   data=dat.tr,method="ML")
    aic.diff<-round(AIC(temp.2)-AIC(temp.1),2)
    temp.3<-data.frame(traitname,x,EStype,c("full","reduced"),
                       c(round(AIC(temp.1),2),round(AIC(temp.2),2)),c("",aic.diff))
    colnames(temp.3)<-colnames(res.aic)
    res.aic<-rbind(res.aic,temp.3)
    #REML for hypothesis testing
    if(aic.diff>2){ #full model if AIC of 2 lower than reduced
      temp.4<-rma.mv(yi,vi,mods=~dat.tr.x+dat.tr.x2,random=list(~1|Species),
                     R=list(Species=phylo.corr),data=dat.tr,method="REML")
      preds<-c(x,paste(x,"2",sep="^"))}else{
        temp.4<-rma.mv(yi,vi,mods=~dat.tr.x,random=list(~1|Species),
                       R=list(Species=phylo.corr),data=dat.tr,method="REML")
        preds<-x}
    temp.6<-data.frame(traitname,x,preds,EStype,round(temp.4$b,5)[-1],
                       round(temp.4$ci.lb,5)[-1],round(temp.4$ci.ub,5)[-1])
    colnames(temp.6)<-colnames(res.tr)
    res.tr<-rbind(res.tr,temp.6)}
  res.tr<-res.tr[-1,];res.aic<-res.aic[-1,]
  res.tr$a<-ifelse(sign(as.numeric(res.tr$CI.l))==sign(as.numeric(res.tr$CI.u)),"*","")
  res.aic$a<-ifelse(as.numeric(res.aic$AIC.diff)>2,"*","")
  res.aic$a[is.na(res.aic$a)]<-""
  if(print.aic==T){res.aic}else{res.tr}}
#Plot trait ES's by isotopes
#Input:
#dat.tr is the output from calc.dat.tr or calc.dat.tr_outer
#traitname is the name of the trait (ex "gill_rakers")
#EStype is "adjacent" or "outer"
#legendPosition says where on the graphs to put genus legend. Usually either "none" or c(0.8,0.8)
#graphs is T or F - to print graphs or table
byIso<-function(dat.tr,traitname,EStype,legendPosition,graphs,phylo){
  
  colnames(dat.tr)[c(3,4,13)]<-paste(colnames(dat.tr)[c(3,4,13)],traitname,sep="_")
  
  #grab the right isotope dataframe for the EStype
  if(EStype=="adjacent SMD"){dat.N15.<-dat.N15}
  if(EStype=="outer SMD"){dat.N15.<-dat.N15.o}
  if(EStype=="lnCVR"){dat.N15.<-dat.N15.cv}
  if(EStype=="adjacent SMD"){dat.C13.<-dat.C13}
  if(EStype=="outer SMD"){dat.C13.<-dat.C13.o}
  if(EStype=="lnCVR"){dat.C13.<-dat.C13.cv}
  
  dat.N15_tr<-merge(dat.N15.,dat.tr,by="group") #put all the data together
  dat.N15_tr<-dat.N15_tr[c(1:15,17,18,27,29)] #get the right columns then name them
  colnames(dat.N15_tr)<-c("group","Lake","yi_N15","vi_isotopes","Country","River_system",
                          "Lat","Lon","logArea","logDepthMax","Lacustrine_morphs","Species","sdi_N15",
                          "Genus","Isotope","yi_tr","vi_tr","sdi_tr","trait.tr")
  dat.N15_tr$Genus[dat.N15_tr$Genus=="Prosopium"]<-"Coregonus"
  
  dat.C13_tr<-merge(dat.C13.,dat.tr,by="group") #same for C13
  dat.C13_tr<-dat.C13_tr[c(1:15,17,18,27,29)]
  colnames(dat.C13_tr)<-c("group","Lake","yi_C13","vi_isotopes","Country","River_system",
                          "Lat","Lon","Area","Depth_max","Lacustrine_morphs","Species","sdi_C13",
                          "Genus","Isotope","yi_tr","vi_tr","sdi_tr","trait.tr")
  dat.C13_tr$Genus[dat.C13_tr$Genus=="Prosopium"]<-"Coregonus"
  
  
  if(graphs==T){
    #set genus colours
    colours<-data.frame(cbind(c(col.cor,col.sal),c("Coregonus","Salvelinus")))
    colours<-colours[colours$X2==levels(as.factor(dat.C13_tr$Genus)),]
    
    #plots! then print plots
    A<-ggplot(data=dat.N15_tr,aes(x=yi_tr,y=yi_N15,colour=Genus))+
      geom_point(size=3)+theme_classic()+
      scale_colour_manual(values=colours$X1)+
      ylab(paste("SMD",EStype,"N15"))+xlab(paste(EStype,traitname))+
      theme(legend.position="none")
    B<-ggplot(data=dat.C13_tr,aes(x=yi_tr,y=yi_C13,colour=Genus))+
      geom_point(size=3)+theme_classic()+
      scale_colour_manual(values=colours$X1)+
      ylab(paste("SMD",EStype,"C13"))+xlab(paste(EStype,traitname))+
      theme(legend.position=legendPosition,
            legend.background=element_rect(linetype=1,size=0.5,colour=1))
    grid.arrange(A,B,ncol=2)}else{
      res.tr<-data.frame("row1","row1","row1","row1","row1","row1")
      
      colnames(res.tr)<-c("isotope","pred","ES_type","estimate","CI.lower","CI.upper")
      #N15
      phylo.corr<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,dat.N15_tr$Species))))
      temp<-rma.mv(yi_N15,vi_isotopes,mods=~yi_tr,random=list(~1|Species),
                   R=list(Species=phylo.corr),data=dat.N15_tr,method="REML")
      temp3<-data.frame("N15",traitname,EStype,round(temp$b,5)[-1],
                        round(temp$ci.lb,5)[-1],round(temp$ci.ub,5)[-1])
      colnames(temp3)<-colnames(res.tr)
      res.tr<-rbind(res.tr,temp3); rm(temp,temp3)
      #C13
      phylo.corr<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,dat.C13_tr$Species))))
      temp<-rma.mv(yi_C13,vi_isotopes,mods=~yi_tr,random=list(~1|Species),
                   R=list(Species=phylo.corr),data=dat.C13_tr,method="REML")
      temp3<-data.frame("C13",traitname,EStype,round(temp$b,5)[-1],
                        round(temp$ci.lb,5)[-1],round(temp$ci.ub,5)[-1])
      colnames(temp3)<-colnames(res.tr)
      res.tr<-rbind(res.tr,temp3)
      res.tr<-res.tr[-1,]
      res.tr$a<-ifelse(sign(as.numeric(res.tr$CI.l))==sign(as.numeric(res.tr$CI.u)),"*","")
      res.tr } }


#calculates difference between two morphs, but with a direction
#Input: tr_cols is the columns numbers of the trait (1) mean (2) n (3) error (4) error type, in that order
#       traitname is a hopefully not useless character string
#       nichename has to be a factor level in the specified nichetype. Like "pelagic" or something
#       nichetype is either "Morph_zone" or "Morph_diet"
#       ES is either "SMD" or "lnCVR"

calc.niche.tr<-function(tr_cols,traitname,nichename,nichetype,ES){
  
  # tr_cols<-106:109
  # traitname<-"body length"
  # nichename<-"littoral/benthic"
  # nichetype<-"Morph_zone"
  # ES<-"lnCVR"
  

  tr<-morphs[c(1:6,tr_cols)]
  colnames(tr)[7:10]<-c("tr.mean","tr.n","tr.e","tr.etype")
  tr$tr.mean<-as.numeric(tr$tr.mean)
  tr$tr.n<-as.numeric(tr$tr.n)
  tr$tr.e<-as.numeric(tr$tr.e)
  tr<-tr[!is.na(tr$tr.n),] # remove if no sample size
  tr<-tr[!is.na(tr$tr.e),] # remove if no error
  #convert Variance, SE and 95CI to SD
  tr$tr.e[tr$tr.etype=="Variance"]<-sqrt(tr$tr.e[tr$tr.etype=="Variance"])
  tr$tr.e[tr$tr.etype=="SE"]<-se_to_sd(tr$tr.e[tr$tr.etype=="SE"],
                                       tr$tr.n[tr$tr.etype=="SE"])
  tr$tr.e[tr$tr.etype=="95CI"]<-CI95_to_sd(tr$tr.e[tr$tr.etype=="95CI"],
                                           tr$tr.n[tr$tr.etype=="95CI"])
  tr$group<-paste(tr$Paper_ID,tr$Lake,sep="_") #make group factor, unique to each lake/study combo
  
  
  tr.z<-tr[tr[,nichetype]!="",] #get dataframe with data for niche type
  tab.z<-data.frame(addmargins(table(tr$group,tr[,nichetype])))
  temp<-tab.z[tab.z$Var2=="Sum"&tab.z$Freq>1,] 
  #keep tab.z groups with 2+ morphs with IDs and data (removing "Sum" group ID and rows)
  tab.z<-tab.z[tab.z$Var1%in%temp$Var1[-length(temp$Var1)]&tab.z$Var2!="Sum",]
  #keep groups with an ID for the niche of interest
  niche<-tr.z[tr.z$group%in%tab.z[tab.z$Var2==nichename&tab.z$Freq>0,]$Var1,]
  niche.tr<-data.frame("row1","row1","row1","row1","row1","row1","row1") #make a new dataframe
  colnames(niche.tr)<-c("group","yi","vi","ID_1","ID_2","n1i","n2i")
  groups<-levels(as.factor(niche$group)) #get a list of groups with relevant data
  
  for(i in groups){
    z<-niche[niche$group==i,c(2,5:9,11)] #select needed columns
    if(nichetype=="Morph_zone"){z<-z[,c(1:2,4:7)]}
    if(nichetype=="Morph_diet"){z<-z[,c(1,3:7)]}
    names(z)<-c("Lake","morphID","mean","n","sd","group") #rename columns
    i1<-which(z$morphID==nichename)
    i2<-which(z$morphID!=nichename)
    
    if(length(i2)==0){print(paste("All morph IDs same for",i))}else{
      if(length(i1)==0){print(paste("No",nichename,"for",i))}else{
        if(length(i1)>0){
          for (j in 1:length(i1)){
            i1a<-i1[j]
            z1<-data.frame(m1i=z$mean[i1a],m2i=z$mean[i2], #make matrix with means, sd's and sample sizes
                           n1i=z$n[i1a],n2i=z$n[i2],sd1i=z$sd[i1a],sd2i=z$sd[i2],
                           group=z$group[i2],ID_1=z$morphID[i1a],ID_2=z$morphID[i2])
            
            #Remove morph if sample size is less than 5
            z1<-z1[z1$n1i>=5,];z1<-z1[z1$n2i>=5,]
            #print error message if no comparison to morph with sample size greater than 5
            if(length(z1$m1i)==0){print(paste("sample size too small for",i))}else{
              if(ES=="SMD"){ #start SMD calc
                z1<-cbind(z1,as.data.frame(escalc(measure="SMD",m1i=z1$m1i,m2i=z1$m2i, #use escalc for SMD and var calculations
                                                  n1i=z1$n1i,n2i=z1$n2i, ## returns the same result as manual eqns below
                                                  sd1i=z1$sd1i,sd2i=z1$sd2i)))
                attr(z1$yi,"measure")<-NULL;attr(z1$yi,"ni")<-NULL #remove attributes
                z1$ID_1<-paste(z1$ID_1,j,sep="")
                niche.tr<-rbind(niche.tr,z1[,c(7,10,11,8,9,3,4)])}
              if(ES=="lnCVR"){ #start lnCVR calc
                yi<-log((z1$sd1i/z1$m1i)/(z1$sd2i/z1$m2i))+ #calc effect size for lnCVR
                  0.5*(1/(z1$n1i-1)-1/(z1$n2i-1))+
                  0.5*(z1$sd2i^2/(z1$n2i*z1$m2i^2)-z1$sd1i^2/(z1$n1i*z1$m1i^2))
                vi<-z1$sd2i^2/(z1$n2i*z1$m2i^2)+ #calc variance for lnCVR
                  z1$sd2i^4/(2*z1$n2i^2*z1$m2i^4)+
                  z1$n2i/(z1$n2i-1)^2+
                  z1$sd1i^2/(z1$n1i*z1$m1i^2)+
                  z1$sd1i^4/(2*z1$n1i^2*z1$m1i^4)+
                  z1$n1i/(z1$n1i-1)^2
                niche.tr.i<-cbind(i,yi,vi,z1$ID_1,z1$ID_2,z1$n1i,z1$n2i)
                colnames(niche.tr.i)<-colnames(niche.tr)
                niche.tr<-rbind(niche.tr,niche.tr.i)}}}
          rm(z1)}}}}
  if(niche.tr$group[1]=="row1"){niche.tr<-niche.tr[-1,]} #Remove useless row
  
  niche.tr<-niche.tr[!is.na(niche.tr$yi),] #only keep rows with data
  
  niche.tr$Lake<-t(data.frame(strsplit(niche.tr$group,"_")))[,2] #add Lake as a variable
  
  temp<-data.frame(table(niche.tr$ID_2,niche.tr$Lake)) #require 5 lakes with data for each morph type
  temp<-temp[temp$Freq>0,] #remove Lake/Morph combos with no data
  temp<-data.frame(table(temp$Var1)) #get counts of lake info by morph
  niche.tr<-niche.tr[niche.tr$ID_2%in%temp[temp$Freq>4,1],] #keep morphs with 5+ observations
  if(nrow(niche.tr)==0){print("STOP: not enough data for this trait")}else{
    if(ES=="SMD"){#Get a combined estimate where there are two "alternate" morphs in a lake - SMD only
      niche.tr$nObs<-"x" #add number of observations column
    
      for(l in levels(as.factor(niche.tr$ID_2))){
        z4<-niche.tr[niche.tr$ID_2==l,]
        for(i in levels(as.factor(z4$group))){
          z2<-z4[z4$group==i,]
          if(nrow(z2)>1){
            z2$yi<-as.numeric(z2$yi);z2$vi<-as.numeric(z2$vi)
            z2$n1i<-as.numeric(z2$n1i);z2$n2i<-as.numeric(z2$n2i)
            #this is the beginning of modified calc.v
            n<-nrow(z2) # number of effect sizes
            V<-matrix(0,nrow=n, ncol=n) # zero matrix
            if(length(levels(as.factor(z2$ID_1)))==1){N<-sum(z2$n1i[1],z2$n2i)} #total sample size
            if(length(levels(as.factor(z2$ID_1)))>1){
              N<-sum(z2[z2$ID_1==levels(as.factor(z2$ID_1))[1],]$n2i)#non-focal morphs sample size
              for(s1 in levels(as.factor(z2$ID_1))){N<-N+z2[z2$ID_1==s1,"n1i"][1]}}
            #calculate variances with total N (if n1 + n2 used in place of N, matches escalc output)
            z2$vi<-1/z2$n1i + 1/z2$n2i+z2$yi^2/(2*N) #from eqn 19.18 Gleser and Olkin 2009
            #calculate covariances where individuals overlap
            #morph of interest is shared individuals
            k<-length(levels(as.factor(z2$ID_1))) #number of ID1 morphs
            m<-length(z2[z2$ID_1==levels(as.factor(z2$ID_1))[1],1]) #number of ID2 morphs
            if(k==1){
              for(j in 1:(n-1)){ # calculate and insert covariances
                v1<-1/z2$n1i[j]+z2$yi[j]*z2$yi[(j+1)]/(2*N) #from eqn 19.19 Gleser and Olkin 2009
                V[j,(j+1)]<-v1;V[(j+1),j]<-v1}}
            #account for shared variance if there are two levels of the focal morph
            if(k>1){ # calculate and insert covariances
              if(m>1){ #if multiple ID2 morphs
                for(j in c(1:(m-1),1:(m-1)+m)){ #get rows that share an ID1
                  v1<-1/z2$n1i[j]+z2$yi[j]*z2$yi[(j+1)]/(2*N) #from eqn 19.19 Gleser and Olkin 2009
                  V[j,(j+1)]<-v1;V[(j+1),j]<-v1}}
              for(o in 1:m){ #get rows that share an ID2
                v1<-1/z2$n2i[o]+z2$yi[o]*z2$yi[(o+m)]/(2*N)
                V[o,(o+m)]<-v1;V[(o+m),o]<-v1}}
            diag(V)<-z2$vi 
            #calculate combined yi
            if(is.numeric(z2$yi)){
              X<-t(t(rep(1,nrow(V))))
              E<-t(t(z2$yi))
              SMD<-solve(t(X)%*%solve(V)%*%X)*(t(X)%*%solve(V)%*%E) #eqn 19.21 from Gleser and Olkin 2009
              #(and eqn 3 from Lajeunesse 2008)
              #calculation has the same result as rma.mv
              vSMD<-solve(t(X)%*%solve(V)%*%X)} #also eqn 19.21 (and paragraph below eqn 3)
            
            z3<-data.frame(cbind(z2$group[1],SMD,vSMD,z2$ID_1[1],z2$ID_2[1],
                                 z2$n1i[1],sum(z2$n2i),z2$Lake[1],n))
            colnames(z3)<-colnames(z2)
            z4<-z4[z4$group!=i,]
            z4<-rbind(z4,z3)}}
        niche.tr<-niche.tr[niche.tr$ID_2!=l,]
        niche.tr<-rbind(niche.tr,z4)}
    niche.tr[niche.tr$nObs=="x","nObs"]<-1}
    
    niche.tr$yi<-as.numeric(niche.tr$yi)
    niche.tr$vi<-as.numeric(niche.tr$vi); niche.tr<-niche.tr[!is.na(niche.tr$vi),]
    
    
    #fixed effects meta-analysis to get one estimate per lake
    tab1<-data.frame(table(niche.tr$Lake, niche.tr$ID_2)) #counts of estimates
    tab1<-tab1[tab1$Freq>1,] #get lakes/ID_2's with repeated estimates
    if(nrow(tab1)>0){ #if there are multiple estimates for a lake
      for(x in 1:nrow(tab1)){ #run through repeated estimate lakes
        dat.sub1<-niche.tr[niche.tr$Lake==tab1$Var1[x]&niche.tr$ID_2==tab1$Var2[x],]
        res.sub1<-rma(dat.sub1$yi,dat.sub1$vi,method="FE") #fixed effects meta-analysis
        #results formatting to add to dat.tr
        Lake<-dat.sub1$Lake[1]
        dattr.sub1<-data.frame(cbind(paste("01E",Lake,sep="_"),res.sub1$beta,
                                     res.sub1$vb, #sqrt(vi)==res.sub1$se
                                     dat.sub1$ID_1[1],dat.sub1$ID_2[1],
                                     sum(as.numeric(dat.sub1$n1i)),
                                     sum(as.numeric(dat.sub1$n2i)),
                                     Lake,dat.sub1$nObs[1]))
        colnames(dattr.sub1)<-colnames(niche.tr); rownames(dattr.sub1)<-NULL
        dattr.sub1$yi<-as.numeric(dattr.sub1$yi); dattr.sub1$vi<-as.numeric(dattr.sub1$vi)
        niche.tr<-niche.tr[!(niche.tr$Lake==tab1$Var1[x]&niche.tr$ID_2==tab1$Var2[x]),] #remove lake from dat.tr
        niche.tr<-rbind(niche.tr,dattr.sub1)}} #add new estimate to dat.tr
    
    #add lake data to dat.gr dataframe
    niche.tr<-merge(niche.tr,df,by="Lake",all.x=T)
    
    niche.tr<-niche.tr[,c("Lake","group","yi","vi","ID_1","ID_2","Country","River_system",
                          "Lat","Lon","logArea","logDepthMax","Lacustrine_morphs")]
    niche.tr<-merge(niche.tr,distinct(niche[,c(3,2)]),by="Lake",all.y=T)
    niche.tr<-niche.tr[!is.na(niche.tr$yi),]
    
    niche.tr$sdi<-sqrt(niche.tr$vi)
    niche.tr$Genus<-t(data.frame(strsplit(niche.tr$Species,"_")))[,1]
    niche.tr$Trait<-traitname
    niche.tr$ID_1<-nichename; niche.tr}}

#dat.tr is the dataframe outputted from calc.dat.tr
#traitname is the name of the trait being added
#meantraits is a preexisting 7 column dataframe with column names/order "mean","CI.l","CI.u","trait","genus","ID_1","ID_2"
#meantraits should have some sort of first row, it will be removed if mean is "row1"
#this function uses rm.mv to do a multivariate meta-analysis of the data for a trait
##with Lake and Genus as random effects if multi=F
##with group, Lake, Genus as random effects if multi=T

add_to_meanNtraits<-function(niche.tr,traitname,IDs,meantraits,multi){

  niche.tr<-niche.tr[!is.na(niche.tr$Lake),]
  for(i in IDs){ #calculate separately for each niche comparison
    
    niche.tr.C<-subset(niche.tr,niche.tr$Genus=="Coregonus"&niche.tr$IDs==i)
    niche.tr.S<-subset(niche.tr,niche.tr$Genus=="Salvelinus"&niche.tr$IDs==i)
    if(multi==T){ #include group as random if multiple traits
      if(dim(niche.tr.C)[1]>2){
        phylo.corr.C<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,niche.tr.C$Species))))
        I2.C<-rma.uni(yi,vi,data=niche.tr.C,method="REML")$I2
        if(phylo.corr.C[,1]!="NaN"){
          rma_tr.C<-rma.mv(yi,vi,random=list(~1|Species,~1|group),R=list(Species=phylo.corr.C),
                           data=subset(niche.tr,niche.tr$Genus=="Coregonus"&niche.tr$IDs==i))}
        if(phylo.corr.C[,1]=="NaN"){
          rma_tr.C<-rma.mv(yi,vi,random=~1|group,
                           data=subset(niche.tr,niche.tr$Genus=="Coregonus"&niche.tr$IDs==i))}}
      if(dim(niche.tr.S)[1]>2){
        phylo.corr.S<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,niche.tr.S$Species))))
        I2.S<-rma.uni(yi,vi,data=niche.tr.S,method="REML")$I2
        if(phylo.corr.S[,1]!="NaN"){
          rma_tr.S<-rma.mv(yi,vi,random=list(~1|Species,~1|group),R=list(Species=phylo.corr.S),
                           data=subset(niche.tr,niche.tr$Genus=="Salvelinus"&niche.tr$IDs==i))}
        if(phylo.corr.S[,1]=="NaN"){
          rma_tr.S<-rma.mv(yi,vi,random=~1|group,
                           data=subset(niche.tr,niche.tr$Genus=="Salvelinus"&niche.tr$IDs==i))}}
    }else{
      #leave out group as random (but still include lake) if multiple traits
      if(dim(niche.tr.C)[1]>2){
        phylo.corr.C<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,niche.tr.C$Species))))
        I2.C<-rma.uni(yi,vi,data=niche.tr.C,method="REML")$I2
        if(phylo.corr.C[,1]!="NaN"){
          rma_tr.C<-rma.mv(yi,vi,random=list(~1|Species),R=list(Species=phylo.corr.C),
                           data=subset(niche.tr,niche.tr$Genus=="Coregonus"&niche.tr$IDs==i))}
        if(phylo.corr.C[,1]=="NaN"){
          rma_tr.C<-rma.mv(yi,vi,
                           data=subset(niche.tr,niche.tr$Genus=="Coregonus"&niche.tr$IDs==i))}}
      if(dim(niche.tr.S)[1]>2){
        phylo.corr.S<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,niche.tr.S$Species))))
        I2.S<-rma.uni(yi,vi,data=niche.tr.S,method="REML")$I2
        if(phylo.corr.S[,1]!="NaN"){
          rma_tr.S<-rma.mv(yi,vi,random=list(~1|Species),R=list(Species=phylo.corr.S),
                           data=subset(niche.tr,niche.tr$Genus=="Salvelinus"&niche.tr$IDs==i))}
        if(phylo.corr.S[,1]=="NaN"){
          rma_tr.S<-rma.mv(yi,vi,
                           data=subset(niche.tr,niche.tr$Genus=="Salvelinus"&niche.tr$IDs==i))}}}
    meantraits<-mutate_all(meantraits,as.character)
    if(dim(subset(niche.tr,niche.tr$Genus=="Coregonus"&niche.tr$IDs==i))[1]>2){
      meantraits<-rbind(meantraits,c(rma_tr.C$beta,rma_tr.C$ci.lb,rma_tr.C$ci.ub,
                                     traitname,"Coregonus",i,I2.C,
                                     nlevels(as.factor(niche.tr.C$Lake))))}
    if(dim(subset(niche.tr,niche.tr$Genus=="Salvelinus"&niche.tr$IDs==i))[1]>2){
      meantraits<-rbind(meantraits,c(rma_tr.S$beta,rma_tr.S$ci.lb,rma_tr.S$ci.ub,
                                     traitname,"Salvelinus",i,I2.S,
                                     nlevels(as.factor(niche.tr.S$Lake))))}}
  if(meantraits$mean[1]=="row1"){meantraits<-meantraits[-1,]}
  meantraits}



#plots mean effect sizes for each trait, or outputs a table
#Variables: tr_cols,traitname,output
#output is "table" or "graph1" or "graph2"
#tr_cols and traitname are same as for calc.dat.tr
#Note that multi is set to F for add_to_meanNtraits
plot.niche.tr<-function(tr_cols,traitname,all.dat,ES){

  #Morph_zone
  pel.tr<-calc.niche.tr(tr_cols,traitname,"pelagic","Morph_zone",ES)
  lit.tr<-calc.niche.tr(tr_cols,traitname,"littoral/benthic","Morph_zone",ES)
  pro.tr<-calc.niche.tr(tr_cols,traitname,"profundal","Morph_zone",ES)
  
  #Morph_diet
  pla.tr<-calc.niche.tr(tr_cols,traitname,"planktivore","Morph_diet",ES)
  ben.tr<-calc.niche.tr(tr_cols,traitname,"benthivore","Morph_diet",ES)
  pis.tr<-calc.niche.tr(tr_cols,traitname,"piscivore","Morph_diet",ES)

  
  datN.all<-rbind(pel.tr,lit.tr,pro.tr,pla.tr,ben.tr,pis.tr)
  
  datN.all<-datN.all[!datN.all$Species%in%c("Prosopium_coulterii","STOP: not enough data for this trait"),]
  
  datN.all[,c("yi","vi")]<-sapply(datN.all[,c("yi","vi")],as.numeric)
  
  datN.all$ID.A<-paste(datN.all$group,datN.all$Genus,datN.all$Trait,
                       datN.all$ID_1,datN.all$ID_2,sep="_")
  datN.all$ID.B<-paste(datN.all$group,datN.all$Genus,datN.all$Trait,
                       datN.all$ID_2,datN.all$ID_1,sep="_")
  
  j<-arrange(datN.all[datN.all$ID.A%in%datN.all$ID.B,],ID_1,ID_2) #select duplicate rows
  for(i1 in levels(as.factor(datN.all$ID.A))){
    if(i1%in%j$ID.B&i1%in%j$ID.A){ #if i1 is in both ID.A and ID.B
      j<-j[i1!=j$ID.B,]}}  #delete the ID.B copy
  datN.all<-rbind(datN.all[!datN.all$ID.A%in%datN.all$ID.B,],j) #combine j and rows that weren't duplicates
  
  datN.all$IDs<-paste(datN.all$ID_1,datN.all$ID_2,sep="-")
  datN.all$ID.A<-NULL;datN.all$ID.B<-NULL
  
  
  #set up pairs of opposite ecotype comparisons; change a1's to corresponding a2's
  a1<-c("piscivore-benthivore","planktivore-benthivore","profundal-littoral/benthic",
        "pelagic-littoral/benthic","planktivore-piscivore")
  a2<-c("benthivore-piscivore","benthivore-planktivore","littoral/benthic-profundal",
        "littoral/benthic-pelagic","piscivore-planktivore")
  
  #switch direction of opposite direction ecotype pairs
  for(i in 1:length(a1)){
    datN.all[datN.all$IDs==a1[i],"yi"]<-
      sapply(datN.all[datN.all$IDs==a1[i],"yi"],byMinusOne)
    datN.all$IDs[datN.all$IDs==a1[i]]<-a2[i]}
  
  
  
  if(all.dat==T){
    datN.all}else{
      #Set up dataframe for mean trait with niches calculation
      datN.all$ID_1<-NULL;datN.all$ID_2<-NULL
      mean_ntraits<-as.data.frame(t(rep("row1",8))) #adjacent dataframe
      colnames(mean_ntraits)<-cbind("mean","CI.l","CI.u","trait","genus",
                                    "IDs","I2","nLakes")
      IDs<-levels(as.factor(datN.all$IDs))
      
      mean_ntraits<-add_to_meanNtraits(datN.all,traitname,IDs,mean_ntraits,F)
      
      mean_ntraits$mean<-as.numeric(mean_ntraits$mean)
      mean_ntraits$CI.l<-as.numeric(mean_ntraits$CI.l);mean_ntraits$CI.u<-as.numeric(mean_ntraits$CI.u)
      mean_ntraits$type<-rep("x",nrow(mean_ntraits))
      mean_ntraits$type[substr(mean_ntraits$IDs,1,3)%in%
                          c("pla","ben","pis")]<-"diet"
      mean_ntraits$type[substr(mean_ntraits$IDs,1,3)%in%
                            c("pel","lit","pro")]<-"habitat"
      mean_ntraits}}




#function to calculation lnCVR for each group, between smallest and largest morph for a given trait
#Input: tr_cols and traitnames are same as for calc.dat.tr
calc.datCV.tr<-function(tr_cols,traitname){
  
  tr<-morphs[c(1:6,tr_cols)]
  colnames(tr)[7:10]<-c("tr.mean","tr.n","tr.e","tr.etype")
  tr$tr.mean<-as.numeric(tr$tr.mean)
  tr$tr.n<-as.numeric(tr$tr.n)
  tr$tr.e<-as.numeric(tr$tr.e)
  tr<-tr[!is.na(tr$tr.n),] # remove if no sample size
  tr<-tr[tr$tr.n>4,] #remove if sample size less than 5
  tr<-tr[!is.na(tr$tr.e),] # remove if no error
  #convert Variance, SE and 95CI to SD
  tr$tr.e[tr$tr.etype=="Variance"]<-sqrt(tr$tr.e[tr$tr.etype=="Variance"])
  tr$tr.e[tr$tr.etype=="SE"]<-se_to_sd(tr$tr.e[tr$tr.etype=="SE"],
                                       tr$tr.n[tr$tr.etype=="SE"])
  tr$tr.e[tr$tr.etype=="95CI"]<-CI95_to_sd(tr$tr.e[tr$tr.etype=="95CI"],
                                           tr$tr.n[tr$tr.etype=="95CI"])
  tr$group<-paste(tr$Paper_ID,tr$Lake,sep="_") #make group factor, unique to each lake/study combo
  groups<-levels(as.factor(tr$group))
  dat.tr<-data.frame("row1","row1","row1");colnames(dat.tr)<-c("group","yi","vi")
  
  for(i in groups){
    z<-tr[tr$group==i,c(2,7:9,11)] #select needed columns
    names(z)<-c("Lake","mean","n","sd","group") #rename columns
    z<-z[order(z$mean),]
    if(nrow(z)<2){print(paste("sample size too small for",i))}else{
      z<-z[c(1,nrow(z)),]
      yi<-log((z$sd[1]/z$mean[1])/(z$sd[2]/z$mean[2]))+
        0.5*(1/(z$n[1]-1)-1/(z$n[2]-1))+
        0.5*(z$sd[2]^2/(z$n[2]*z$mean[2]^2)-z$sd[1]^2/(z$n[1]*z$mean[1]^2))
      vi<-z$sd[2]^2/(z$n[2]*z$mean[2]^2)+
        z$sd[2]^4/(2*z$n[2]^2*z$mean[2]^4)+
        z$n[2]/(z$n[2]-1)^2+
        z$sd[1]^2/(z$n[1]*z$mean[1]^2)+
        z$sd[1]^4/(2*z$n[1]^2*z$mean[1]^4)+
        z$n[1]/(z$n[1]-1)^2
      x<-cbind(i,yi,vi)
      colnames(x)<-colnames(dat.tr)
      dat.tr<-rbind(dat.tr,x)}}
  dat.tr<-dat.tr[-1,]
  dat.tr$group<-as.character(dat.tr$group)
  dat.tr$yi<-as.numeric(as.character(dat.tr$yi))
  dat.tr$vi<-as.numeric(as.character(dat.tr$vi))
  dat.tr$Lake<-t(data.frame(strsplit(dat.tr$group,"_")))[,2]
  
  
  #fixed effects meta-analysis to get one estimate per lake
  tab1<-data.frame(table(dat.tr$Lake)) #counts of estimates
  tab1<-tab1[tab1$Freq>1,] #get lakes with repeated estimates
  if(nrow(tab1)>0){ #if there are multiple estimates for a lake
    for(x in 1:nrow(tab1)){ #run through repeated estimate lakes
      dat.sub1<-dat.tr[dat.tr$Lake==tab1$Var1[x],]
      res.sub1<-rma(dat.sub1$yi,dat.sub1$vi,method="FE") #fixed effects meta-analysis
      #results formatting to add to dat.tr
      Lake<-dat.sub1$Lake[1]
      dattr.sub1<-data.frame(cbind(paste("01E",Lake,sep="_"),res.sub1$beta,
                                   res.sub1$vb, #sqrt(vi)==res.sub1$se
                                   Lake))
      colnames(dattr.sub1)<-colnames(dat.tr); rownames(dattr.sub1)<-NULL
      dattr.sub1$yi<-as.numeric(dattr.sub1$yi); dattr.sub1$vi<-as.numeric(dattr.sub1$vi)
      dat.tr<-dat.tr[dat.tr$Lake!=tab1$Var1[x],] #remove lake from dat.tr
      dat.tr<-rbind(dat.tr,dattr.sub1)}} #add new estimate to dat.tr
  dat.tr[dat.tr$Lake=="Kronotskoe",]
  #add lake data to dat.gr dataframe
  dat.tr<-merge(dat.tr,df,by="Lake",all.x=T)
  dat.tr<-dat.tr[,c("Lake","group","yi","vi","Country","River_system",
                    "Lat","Lon","logArea","logDepthMax","Lacustrine_morphs",
                    "Temp_mean","Temp_sea","Prec_mean","Prec_sea",
                    "Lat.2","logArea.2","logDepthMax.2",
                    "Temp.2","Prec.2","Continent",
                    "PC1_climate","PC2_size","PC3_depthRatio","PC4_Latitude")]
  dat.tr<-merge(dat.tr,distinct(tr[,c(3,2)]),by="Lake",all.y=T)
  if(nrow(dat.tr)==length(unique(tr$Lake))){print("merging probably worked okay")}else{
    print("something went wrong, check merging")}
  dat.tr$sdi<-sqrt(dat.tr$vi)
  dat.tr$Genus<-t(data.frame(strsplit(dat.tr$Species,"_")))[,1]
  dat.tr$Trait<-traitname
  dat.tr<-subset(dat.tr,!is.na(dat.tr$vi)); dat.tr}

#niche.tr is the dataframe outputted from calc.niche.tr
#traitname is the name of the trait being added
#nichename is the niche of interest
#phylo is a Newick format phylogeny
#meantraits is a preexisting 7 column dataframe with column names/order "mean","CI.l","CI.u","trait","genus","ID_1","ID_2"
#meantraits should have some sort of first row, it will be removed if mean is "row1"
#this function uses rm.mv to do a multivariate meta-analysis of the data for a trait

add_to_meanNCVtraits<-function(niche.tr,traitname,nichename,meantraits,phylo){
  niche.tr<-niche.tr[!is.na(niche.tr$Lake),]
  niche.tr$Genus[niche.tr$Genus=="Prosopium"]<-"Coregonus"
  for(i in levels(as.factor(niche.tr$ID_2))){ #calculate separately for each niche category
    #leave out group as random (but still include lake) if multiple traits
    #Coregonus
    niche.tr.C<-subset(niche.tr,niche.tr$Genus=="Coregonus"&niche.tr$ID_2==i)
    if(length(levels(as.factor(niche.tr.C$Species)))>1){
      phylo.corr.C<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,niche.tr.C$Species))))
      if(dim(niche.tr.C)[1]>2){
        rma_tr.C<-rma.mv(yi,vi,random=list(~1|Species,~1|Lake,~1|group),
                         R=list(Species=phylo.corr.C),data=niche.tr.C)}}
    if(length(levels(as.factor(niche.tr.C$Species)))==1){
      if(dim(niche.tr.C)[1]>2){rma_tr.C<-rma.mv(yi,vi,random=list(~1|Lake,~1|group),data=niche.tr.C)}}
    #Salvelinus
    niche.tr.S<-subset(niche.tr,niche.tr$Genus=="Salvelinus"&niche.tr$ID_2==i)
    if(length(levels(as.factor(niche.tr.S$Species)))>1){
      phylo.corr.S<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,niche.tr.S$Species))))
      if(dim(niche.tr.S)[1]>2){
        rma_tr.S<-rma.mv(yi,vi,random=list(~1|Species,~1|Lake,~1|group),
                         R=list(Species=phylo.corr.S),data=niche.tr.S)}}
    if(length(levels(as.factor(niche.tr.S$Species)))==1){
      if(dim(niche.tr.S)[1]>2){rma_tr.S<-rma.mv(yi,vi,random=list(~1|Lake,~1|group),data=niche.tr.S)}}
    meantraits<-mutate_all(meantraits,as.character)
    if(dim(niche.tr.C)[1]>2){
      meantraits<-rbind(meantraits,c(rma_tr.C$beta,rma_tr.C$ci.lb,rma_tr.C$ci.ub,
                                     traitname,"Coregonus",nichename,i))}
    if(dim(niche.tr.S)[1]>2){
      meantraits<-rbind(meantraits,c(rma_tr.S$beta,rma_tr.S$ci.lb,rma_tr.S$ci.ub,
                                     traitname,"Salvelinus",nichename,i))}}
  if(meantraits$mean[1]=="row1"){meantraits<-meantraits[-1,]}
  meantraits}


#Plots or produces a table for comparison of lnCVR between niche categories
plot.nicheCV.tr<-function(tr_cols,traitname,all.dat){

  #Morph_zone
  pel.tr<-calc.niche.tr(tr_cols,traitname,"pelagic","Morph_zone","lnCVR")
  lit.tr<-calc.niche.tr(tr_cols,traitname,"littoral/benthic","Morph_zone","lnCVR")
  pro.tr<-calc.niche.tr(tr_cols,traitname,"profundal","Morph_zone","lnCVR")
  
  #Morph_diet
  pla.tr<-calc.nicheCV.tr(tr_cols,traitname,"planktivore","Morph_diet","lnCVR")
  ben.tr<-calc.nicheCV.tr(tr_cols,traitname,"benthivore","Morph_diet","lnCVR")
  pis.tr<-calc.nicheCV.tr(tr_cols,traitname,"piscivore","Morph_diet","lnCVR")
  
  if(all.dat==T){
    mean_ntraits<-rbind(pel.tr,lit.tr,pro.tr,pla.tr,ben.tr,pis.tr)
    #mean_ntraits$genus<-t(data.frame(strsplit(mean_ntraits$Species,"_")))[,1]
  }else{
  
  #Set up dataframe for mean trait with niches calculation
  mean_ntraits<-as.data.frame(t(rep("row1",7))) #adjacent dataframe
  colnames(mean_ntraits)<-cbind("mean","CI.l","CI.u","trait","genus","ID_1","ID_2")
  
  if(is.data.frame(pel.tr)){mean_ntraits<-add_to_meanNCVtraits(pel.tr,traitname,"pelagic",mean_ntraits,phylo)}
  if(is.data.frame(lit.tr)){mean_ntraits<-add_to_meanNCVtraits(lit.tr,traitname,"littoral/benthic",mean_ntraits,phylo)}
  if(is.data.frame(pro.tr)){mean_ntraits<-add_to_meanNCVtraits(pro.tr,traitname,"profundal",mean_ntraits,phylo)}
  
  if(is.data.frame(pla.tr)){mean_ntraits<-add_to_meanNCVtraits(pla.tr,traitname,"planktivore",mean_ntraits,phylo)}
  if(is.data.frame(ben.tr)){mean_ntraits<-add_to_meanNCVtraits(ben.tr,traitname,"benthivore",mean_ntraits,phylo)}
  if(is.data.frame(pis.tr)){mean_ntraits<-add_to_meanNCVtraits(pis.tr,traitname,"piscivore",mean_ntraits,phylo)}
  
  mean_ntraits$mean<-as.numeric(mean_ntraits$mean)
  mean_ntraits$CI.l<-as.numeric(mean_ntraits$CI.l);mean_ntraits$CI.u<-as.numeric(mean_ntraits$CI.u)
  mean_ntraits<-mean_ntraits[!duplicated(abs(mean_ntraits$mean)),]}
  mean_ntraits$type<-rep("x",nrow(mean_ntraits))
  mean_ntraits$type[mean_ntraits$ID_1%in%
                      c("planktivore","benthivore","piscivore")]<-"diet"
  mean_ntraits$type[mean_ntraits$ID_1%in%
                      c("pelagic","littoral/benthic","profundal")]<-"habitat"
  
  # mean_ntraits$ID.A<-paste(mean_ntraits$genus,paste(mean_ntraits$ID_1,
  #                                                   mean_ntraits$ID_2,sep="_"),sep="_")
  # mean_ntraits$ID.B<-paste(mean_ntraits$genus,paste(mean_ntraits$ID_2,
  #                                                   mean_ntraits$ID_1,sep="_"),sep="_")
  # x<-mean_ntraits[mean_ntraits$ID.A%in%mean_ntraits$ID.B,]
  # y<-x
  # for(i in 1:(nrow(x)/2)){
  #   x<-x[-which(x$ID.A[i]==x$ID.B),]}
  # mean_ntraits<-mean_ntraits[!mean_ntraits$ID.A%in%setdiff(y$ID.A,x$ID.A),]
  # mean_ntraits$IDs<-paste(mean_ntraits$ID_1,mean_ntraits$ID_2,sep="-")
  mean_ntraits}


#Calculates observed chi square (lazily) and simulates null distribution for comparison
#Input: Matrix or dataframe with numeric values
#Output: observed chi square, BCa confidence interval, p value based on simulations, histogram simulated chi stats
simulate.chisq<-function(mat.x){
  
  a<-chisq.test(mat.x) #calculate observed stat
  #Randomly generate some tables using Patefield 1981 algorithm
  null.x<-r2dtable(10000,rowSums(mat.x),colSums(mat.x))
  
  z<-chisq.test(null.x[[1]]) #calculate first chisq
  y<-z$statistic #save to object y
  for(x in 2:10000){ #for each simulated matrix
    z<-chisq.test(null.x[[x]]) #calculate chisq
    y<-c(y,z$statistic)} #add chisq to object y
  b<-bca(y,conf.level=0.95) #function from library coxed
  b<-paste(b[1],b[2],sep=",")
  d<-y[y>a$statistic]
  p.value<-1-length(d)/length(y)
  print(paste("Observed stat:",a$statistic)) #report test stat
  print(paste("Chi/df:",a$statistic/a$parameter))
  print(paste("df:",a$parameter))
  #print(paste("BCa Confidence interval:",b)) #report CI
  print(paste("simulated p-value:",p.value))
  hist(y,main="Simulated chi stats")} #histogram of simulated stats

#format data to run k-sample Anderson-Darling test
#morphs.tr is morphs dataframe formatted for trait of interest
#trait_col is the name of the column of trait of interest means
#genus.ID is "Sal" for Salvelinus or "Cor" for Coregonus
#if adTest=T the performs A-D test, otherwise outputs formatted dataframe for plotting
#before running tests, calculates mean values for each ecotype from each lake with multiple samples
phenotype_repeat<-function(morphs.tr,trait_col,n_col,genus.ID,adTest){
  
  # #TEMP TEST DATA
  # morphs.tr<-morphs.fl
  # trait_col<-"Total_length"
  # n_col<-"length_n"
  # genus.ID<-"Sal"
  # adTest<-T
  
  

  morphs.tr$ID<-paste(morphs.tr$Paper_ID,morphs.tr$Lake,sep="_") #make dataframe for trait analysis
  morphs.tr<-morphs.tr[substr(morphs.tr$Species,1,3)==genus.ID,] #keep only one genus
  #start IDing lakes and selecting replicate
  morphs.tr.ID<-data.frame(table(morphs.tr$ID)) #get paper IDS
  morphs.tr.ID$Lake<-t(data.frame(strsplit(as.character(morphs.tr.ID$Var1),"_")))[,2] #extract lake name from paper ID
  morphs.tr.ID<-merge(morphs.tr.ID,df,by="Lake")[,c("Var1","Freq","Lake","Lacustrine_morphs")] #get # of morphs per lake
  morphs.tr.ID$Var1<-as.character(morphs.tr.ID$Var1)
  morphs.tr.ID<-morphs.tr.ID[morphs.tr.ID$Freq==morphs.tr.ID$Lacustrine_morphs,]#only keep cases where all morphs have a mean recorded
  lakes.rep<-as.character(data.frame(table(morphs.tr.ID$Lake)[table(morphs.tr.ID$Lake)>1])$Var1) #find lakes with repeated measures
  morphs.tr<-morphs.tr[morphs.tr$ID%in%morphs.tr.ID$Var1,] #only keep selected IDs
  
  morphs.tr<-merge(morphs.tr,df,by="Lake") #add info about lakes
  #only keep data from papers where all morphs have been identified
  temp2<-morphs.tr[substr(morphs.tr$Species,1,3)==genus.ID, #only keep one genus
                   c("ID","Lake",trait_col,"Species","Lacustrine_morphs","Morph_zone","Morph_diet")]
  
  colnames(temp2)[3]<-"trait"
  temp2<-arrange(temp2,Lake,ID,trait)
  # for(i in 1:length(lakes.rep)){ #pick one randomly from lakes with repeated measures
  #   Var1s<-morphs.tr.ID[morphs.tr.ID$Lake==lakes.rep[i],"Var1"]
  #   morphs.tr.ID<-morphs.tr.ID[!morphs.tr.ID$Var1%in%Var1s[Var1s!=sample(Var1s,1)],]
  # }
  if(length(lakes.rep)>0){
    for(i in 1:length(lakes.rep)){ 
      Var1s<-morphs.tr.ID[morphs.tr.ID$Lake==lakes.rep[i],"Var1"]
      Var.dat<-morphs.tr[morphs.tr$ID%in%Var1s,c("ID",trait_col,n_col)]
      colnames(Var.dat)<-c("ID","trait","tr_n")
      Var.dat<-arrange(Var.dat,ID,trait)
      n<-length(Var.dat$ID)/nlevels(as.factor(Var.dat$ID))
      k<-vector()
      for(j in 1:n){ #make vector with means of each ecotype
        k.i<-seq(j,nrow(Var.dat),n)
        k.mean<-sum(Var.dat$trait[k.i]*Var.dat$tr_n[k.i])/sum(Var.dat$tr_n[k.i])
        k<-c(k,k.mean)}
      temp2<-temp2[temp2$Lake!=lakes.rep[1],]
      temp2<-rbind(temp2,data.frame(ID=paste("01E",lakes.rep[i],sep="_"),
                                    Lake=lakes.rep[i],trait=k,
                                    Species=temp2[temp2$Lake==lakes.rep[i],"Species"][1],
                 Lacustrine_morphs=temp2[temp2$Lake==lakes.rep[i],"Lacustrine_morphs"][1],
                 Morph_zone=temp2[temp2$Lake==lakes.rep[i],"Morph_zone"][1:n],
                 Morph_diet=temp2[temp2$Lake==lakes.rep[i],"Morph_diet"][1:n]))}}
  temp3<-temp2[grep("01E",temp2$ID),]
  temp2<-temp2[!temp2$Lake%in%lakes.rep,]
  temp2<-rbind(temp2,temp3); rm(temp3)
  temp2$ID2<-paste(temp2$ID,1:nrow(temp2),sep="_")
  
  for(i in unique(temp2$ID)){
    df.i<-temp2[temp2$ID==i,]
    df.i$trait-mean(df.i$trait)
  }
  
 
  
  temp1<-as.data.frame(pivot_wider(temp2,names_from=Lake,values_from=trait)) #means in columns with each lake as a column name
  temp1
  if(adTest==T){x<-ad.test(temp1[,7:ncol(temp1)],method="simulated")
                print(x)
                print(1-x$ad[1,3])}else{ #run ad test on numeric columns
    temp2 } }


byMinusOne<-function(x){x*(-1)} #multiply by -1


#runs models to estimate genus, trait, and interaction effects for pairs of ecotypes
#datN.all is output from plot.niche.tr or plot.nicheCV.tr with all.dat=T
calc.datN.all.out<-function(datN.all,ES){
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ES<-"SMD" #temp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  #get vector of ecotype comparisons
  IDs<-levels(as.factor(datN.all$IDs))
  
  
  #make output dataframe
  datN.all.out<-data.frame(ecotypes="row1",mod="row1",estimate="row1",Z.value="row1",
                           p.value="row1",I2="row1",row.names=NULL)
  
  
  datN.ID.plots<-list()
  
  j<-0 #index value for adding plots to list
  
  for(i in 1:length(IDs)){
    
    datN.ID<-datN.all[datN.all$IDs==IDs[i],]
    if(nlevels(as.factor(datN.ID$Genus))>1&nlevels(as.factor(datN.ID$Trait))>1){
      phylo.corr.ID<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,datN.ID$Species))))
      x1<-rma.mv(yi,vi,mods=~Trait*Genus-1,random=list(~1|Species,~1|group),
                 R=list(Species=phylo.corr.ID),data=datN.ID)
      y1<-rma.uni(yi,vi,data=datN.ID)
      datN.all.out<-rbind(datN.all.out,data.frame(ecotypes=IDs[i],mod=colnames(x1$vb),
                                                  estimate=x1$beta,Z.value=x1$zval,
                                                  p.value=x1$pval,I2=y1$I2,
                                                  row.names=NULL))}
    if(nlevels(as.factor(datN.ID$Genus))>1&nlevels(as.factor(datN.ID$Trait))<=1){
      phylo.corr.ID<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,datN.ID$Species))))
      x1<-rma.mv(yi,vi,mods=~Genus-1,random=list(~1|Species,~1|group),
                 R=list(Species=phylo.corr.ID),data=datN.ID)
      y1<-rma.uni(yi,vi,data=datN.ID)
      datN.all.out<-rbind(datN.all.out,data.frame(ecotypes=IDs[i],mod=colnames(x1$vb),
                                                  estimate=x1$beta,Z.value=x1$zval,
                                                  p.value=x1$pval,I2=y1$I2,
                                                  row.names=NULL))}
    if(nlevels(as.factor(datN.ID$Genus))<=1&nlevels(as.factor(datN.ID$Trait))>1){
      phylo.corr.ID<-vcv(compute.brlen(drop.tip(phylo,setdiff(phylo$tip.label,datN.ID$Species))))
      x1<-rma.mv(yi,vi,mods=~Trait-1,random=list(~1|Species,~1|group),
                 R=list(Species=phylo.corr.ID),data=datN.ID)
      y1<-rma.uni(yi,vi,data=datN.ID)
      datN.all.out<-rbind(datN.all.out,data.frame(ecotypes=IDs[i],mod=colnames(x1$vb),
                                                  estimate=x1$beta,Z.value=x1$zval,
                                                  p.value=x1$pval,I2=y1$I2,
                                                  row.names=NULL))}
    
    if((nlevels(as.factor(datN.ID$Genus))+nlevels(as.factor(datN.ID$Trait)))>2){
    
    datN.ID.plot<-datN.ID
    datN.ID.plot$sdi<-as.numeric(datN.ID.plot$sdi)
    datN.ID.plot$CI_l<-datN.ID.plot$yi-1.96*sqrt(datN.ID.plot$vi)
    datN.ID.plot$CI_u<-datN.ID.plot$yi+1.96*sqrt(datN.ID.plot$vi)
    datN.ID.plot<-datN.ID.plot[order(datN.ID.plot$Lake),]
    
    genus.rows<-grep("Genus",rownames(x1$b))
    
    if(length(genus.rows)<length(x1$b)){
    
      
      j<-j+1
    if(length(genus.rows)>0){
      temp1<-data.frame(paste("RE Model:",gsub("Trait","",rownames(x1$b)[-genus.rows])),
                        NA,x1$b[-genus.rows],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                        "both",gsub("Trait","",rownames(x1$b)[-genus.rows]),datN.ID.plot$IDs[1],
                        x1$ci.lb[-genus.rows],x1$ci.ub[-genus.rows])}else{
                          temp1<-data.frame(paste("RE Model:",gsub("Trait","",rownames(x1$b))),
                                            NA,x1$b,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                                            "both",gsub("Trait","",rownames(x1$b)),datN.ID.plot$IDs[1],
                                            x1$ci.lb,x1$ci.ub)
                        }
      colnames(temp1)<-colnames(datN.ID.plot)
      datN.ID.plot<-rbind(datN.ID.plot,temp1)
      
      
      
      datN.ID.plot$Lake<-factor(datN.ID.plot$Lake,levels=rev(unique(datN.ID.plot$Lake)))
      
      
      datN.ID.plots[[j]]<-ggplot(data=datN.ID.plot,aes(x=yi,y=Lake,xmin=CI_l,xmax=CI_u,colour=Genus))+
        geom_vline(xintercept=0,lty=2,colour="grey20")+
        geom_point(aes(size=Genus))+geom_errorbar(width=0.2)+
        facet_grid(Trait~.,scales= "free", space="free")+
        theme_classic()+xlab(paste(ES,datN.ID.plot$IDs[1]))+
        scale_size_manual(values=c(2,1,1))+
        scale_color_manual(values=c("grey20",col.cor,col.sal))}}}
  
  #PROBLEM HERE: means fitted based on old modelling approach
  ##redo with new combined modelling
  
  
  
  
  
  
  
  
  
  datN.all.out<-datN.all.out[datN.all.out$ecotypes!="row1",]
  datN.all.out[,c("estimate","Z.value","p.value","I2")]<-
    sapply(datN.all.out[,c("estimate","Z.value","p.value","I2")],as.numeric)
  
  
  datN.all.out$a<-rep("",nrow(datN.all.out))
  datN.all.out$a[datN.all.out$p.value<0.05]<-"*"
  
  datN.all.out[,c("estimate","Z.value","p.value")]<-round(datN.all.out[,c("estimate","Z.value","p.value")],digits=3)
  datN.all.out$p.value[datN.all.out$p.value==0]<-"<0.001"
  datN.ID.plots[[(j+1)]]<-datN.all.out
  datN.ID.plots}


#Uses output from "dredge" function to calculate AICc weights for a set of 
##predictor variables
#dLm = dredge output
#genus = genus name (Coregonus or Salvelinus)
#vars_list = vector of names of variables that were included as predictors in the model
# calc_AICc_weights<-function(dLm,genus,vars_list){
#   
#   tr_weights<-data.frame()
#   
#   for(tr in vars_list){
#     tr_weight<-sum(dLm$weight[!is.na(dLm[,tr])])
#     tr_weights<-rbind(tr_weights,c(tr,tr_weight,genus))}
#   colnames(tr_weights)<-c("abiotic_var","AICc_weight","Genus")
#   tr_weights$AICc_weight<-as.numeric(tr_weights$AICc_weight)
#   tr_weights}


#Uses output from "glmulti" function to calculate AICc weights for a set of 
##predictor variables
#lm.tr = glmulti output
#Genus = genus name (Coregonus or Salvelinus)
#tr_name = traitname (ex. body_length)
calc_AICc_weights<-function(lm.tr,Genus.i,tr_name){
  
  #lm.tr<-lm.cor.fl
  #Genus.i<-"Coregonus"
  #tr_name<-"body_length"
  
  
  lm.tr.w<-weightable(lm.tr) #save weights for each model
  #delete everything except model terms, split terms using + as pattern
  #saves as a matrix with 4 columns - convert to dataframe and add weights col
  lm.tr.p<-data.frame(str_split_fixed( 
    substr(lm.tr.w$model,10,nchar(lm.tr.w$model))," \\+ ",n=4),
    weights=lm.tr.w$weights)
  
  out.fl<-data.frame() #make output dataframe
  for(pVar in lm.tr@call$xr){ #extract predictors from model and cycle through
    #function uses grep, so check that no var names are subsets of other var names
    if(length(lm.tr@call$xr[grepl(pVar,lm.tr@call$xr)])>1){
      print("PROBLEM: variable name subset of another variable name")}
    df.temp<-data.frame(genus=Genus.i,trait=tr_name, #save genus and trait name
                        pred_var=pVar, #save abiotic var name
                        #select all rows containing a trait name and sum those weights
                        AICc_weight=lm.tr.p%>%filter_all(any_vars(.==pVar))%>%
                          summarise(weight_sum=sum(weights)))
    out.fl<-rbind(out.fl,df.temp) #add results to dataframe
  }
  out.fl}


fit_ordinal_models<-function(df.genus){
  
  brm1<-brm(data = df.genus,family = sratio("cloglog"),
            ordinal_morphs ~ PC4_Latitude)
  brm2<-brm(data = df.genus,family = sratio("cloglog"),
            ordinal_morphs ~ PC2_size + PC4_Latitude)
  brm3<-brm(data = df.genus,family = sratio("cloglog"),
            ordinal_morphs ~ PC1_climate + PC4_Latitude)
  brm4<-brm(data = df.genus,family = sratio("cloglog"),
            ordinal_morphs ~ PC3_depthRatio + PC4_Latitude)
  brm5<-brm(data = df.genus,family = sratio("cloglog"),
            ordinal_morphs ~ PC2_size + PC3_depthRatio + PC4_Latitude)
  brm6<-brm(data = df.genus,family = sratio("cloglog"),
            ordinal_morphs ~ PC1_climate + PC2_size + PC4_Latitude)
  brm7<-brm(data = df.genus,family = sratio("cloglog"),
            ordinal_morphs ~ PC1_climate + PC3_depthRatio + PC4_Latitude)
  brm8<-brm(data = df.genus,family = sratio("cloglog"),
            ordinal_morphs ~ PC1_climate + PC2_size + PC3_depthRatio + PC4_Latitude)
  brm9<-brm(data = df.genus,family = sratio("cloglog"),
            ordinal_morphs ~ 1)
  brm10<-brm(data = df.genus,family = sratio("cloglog"),
             ordinal_morphs ~ PC1_climate)
  brm11<-brm(data = df.genus,family = sratio("cloglog"),
             ordinal_morphs ~ PC2_size)
  brm12<-brm(data = df.genus,family = sratio("cloglog"),
             ordinal_morphs ~ PC3_depthRatio)
  brm13<-brm(data = df.genus,family = sratio("cloglog"),
             ordinal_morphs ~ PC1_climate + PC2_size)
  brm14<-brm(data = df.genus,family = sratio("cloglog"),
             ordinal_morphs ~ PC1_climate + PC3_depthRatio)
  brm15<-brm(data = df.genus,family = sratio("cloglog"),
             ordinal_morphs ~ PC2_size + PC3_depthRatio)
  brm16<-brm(data = df.genus,family = sratio("cloglog"),
             ordinal_morphs ~ PC1_climate + PC2_size + PC3_depthRatio)
  
  
  brmList<-list(brm1,brm2,brm3,brm4,brm5,brm6,brm7,brm8,brm9,brm10,
                brm11,brm12,brm13,brm14,brm15,brm16)
  
  brmList
  
  
  
}

# fit_ordinal_models<-function(df.genus){
#   
#   brm1<-brm(data = df.genus,family = cumulative(probit),
#                ordinal_morphs ~ PC4_Latitude)
#   brm2<-brm(data = df.genus,family = cumulative(probit),
#             ordinal_morphs ~ PC2_size + PC4_Latitude)
#   brm3<-brm(data = df.genus,family = cumulative(probit),
#             ordinal_morphs ~ PC1_climate + PC4_Latitude)
#   brm4<-brm(data = df.genus,family = cumulative(probit),
#             ordinal_morphs ~ PC3_depthRatio + PC4_Latitude)
#   brm5<-brm(data = df.genus,family = cumulative(probit),
#             ordinal_morphs ~ PC2_size + PC3_depthRatio + PC4_Latitude)
#   brm6<-brm(data = df.genus,family = cumulative(probit),
#             ordinal_morphs ~ PC1_climate + PC2_size + PC4_Latitude)
#   brm7<-brm(data = df.genus,family = cumulative(probit),
#             ordinal_morphs ~ PC1_climate + PC3_depthRatio + PC4_Latitude)
#   brm8<-brm(data = df.genus,family = cumulative(probit),
#             ordinal_morphs ~ PC1_climate + PC2_size + PC3_depthRatio + PC4_Latitude)
#   brm9<-brm(data = df.genus,family = cumulative(probit),
#             ordinal_morphs ~ 1)
#   brm10<-brm(data = df.genus,family = cumulative(probit),
#             ordinal_morphs ~ PC1_climate)
#   brm11<-brm(data = df.genus,family = cumulative(probit),
#              ordinal_morphs ~ PC2_size)
#   brm12<-brm(data = df.genus,family = cumulative(probit),
#              ordinal_morphs ~ PC3_depthRatio)
#   brm13<-brm(data = df.genus,family = cumulative(probit),
#              ordinal_morphs ~ PC1_climate + PC2_size)
#   brm14<-brm(data = df.genus,family = cumulative(probit),
#              ordinal_morphs ~ PC1_climate + PC3_depthRatio)
#   brm15<-brm(data = df.genus,family = cumulative(probit),
#              ordinal_morphs ~ PC2_size + PC3_depthRatio)
#   brm16<-brm(data = df.genus,family = cumulative(probit),
#              ordinal_morphs ~ PC1_climate + PC2_size + PC3_depthRatio)
#   
#   
#   brmList<-list(brm1,brm2,brm3,brm4,brm5,brm6,brm7,brm8,brm9,brm10,
#                 brm11,brm12,brm13,brm14,brm15,brm16)
#   
#   brmList
#   
#   
#   
# }
