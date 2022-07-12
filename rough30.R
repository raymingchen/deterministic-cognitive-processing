######## rough sets calculations
rm(list=ls());
library("igraph");
## convert a data vector (an information system) into a partition (an equivalence relation)
EqClass=function(dataColVec)  ## dataColVec:a column vector of data, for example, a property valued over a set of objects
{
rows=length(dataColVec);
uqVec=unique(dataColVec); ## the unique elements in dataColVec
cols=length(uqVec);
eqClass=matrix(0,nrow=rows,ncol=cols);
  for(i in 1:cols)
  {
  posInd=which(dataColVec==uqVec[i]);
  eqClass[posInd,i]=1;  
  }
return(eqClass);
}

### convert a target set into a target column vector
TargetColVec=function(targetSet,lenDom) ##lenDom:length of domain; targetSet is an input in the form of a vector
{
targetColVec=rep(0,lenDom);
targetColVec[targetSet]=1;
return(targetColVec);
}

## compute lower and upper vector (approximatioin), given a target set
LowUpperVecs=function(targetSet,eqClass)   ##eqClass:equivalence class, an |object|x|domain of a feature| matrix
{
lenDom=dim(eqClass)[1];
targetColVec=TargetColVec(targetSet,lenDom);
nOBJ=length(targetColVec);
rowSum=apply(eqClass,2,sum);
IntSec=targetColVec*eqClass;
rowSum_IntSec=apply(IntSec,2,sum);

lowerVecs=rep(0,nOBJ);
upPosInd=which(rowSum_IntSec!=0);
upperVecs=apply(as.matrix(eqClass[,upPosInd]),1,sum);
lowerPosInd=which(rowSum_IntSec==rowSum);
lowerVecs=apply(as.matrix(eqClass[,lowerPosInd]),1,sum);
lowerUpperVecs=matrix(nrow=nOBJ,ncol=2);
lowerUpperVecs[,1]=lowerVecs;
lowerUpperVecs[,2]=upperVecs;
return(lowerUpperVecs)
}


## features representations of a target set
FeaRep=function(targetSet,dataMatrix)   ## |object|-by-|features| matrix
{
lenDom=dim(dataMatrix)[1];
num_features=dim(dataMatrix)[2];
feaRep=as.list(numeric(1*num_features));dim(feaRep)=c(1,num_features);
   for(i in 1:num_features)
   {
   dataColVec=dataMatrix[,i];
   eqClass=EqClass(dataColVec);
   lowUpperVecs=LowUpperVecs(targetSet,eqClass);
   feaRep[[1,i]]=lowUpperVecs;
   }
return(feaRep);
}


### compute distance between feature (target sets) representations
Dis_FeaReps=function(feaRep1,feaRep2) ## two |object|-by-|features| lists
{
num_fea=length(feaRep1);
dis=rep(0,num_fea);
for(i in 1:num_fea)
{
mat1=feaRep1[[i]];
mat2=feaRep2[[i]];
difMat=mat1-mat2;
dis[i]=sum(abs(mat1-mat2));
}
dis_FeaReps=sum(dis);
return(dis_FeaReps);
}

###position for the directed minimal pairs
MinPosInd=function(disMatrix)
{
n=dim(disMatrix)[1];
min=apply(disMatrix,1, FUN=function(x){min(x[x>0])});
minPosInd=diag(n);
  for(i in 1:n)
  {
  veci=disMatrix[i,];
  ind=which(veci==min[i]);
  minPosInd[i,ind]=1;  
  }
return(minPosInd);
}


###check symmetric elements (bi-directed minimal pairs)
BiDirMinPairs=function(minPosInd)
{
biDirMinPairs=minPosInd>0 & minPosInd==t(minPosInd);
biDirMinPairs=biDirMinPairs*1;  ##convert true/false to 0/1
return(biDirMinPairs);
}


### targets congregated
TarCong=function(biDirMinPairs)   ##bi-directed minimal pairs
{
graph=graph_from_adjacency_matrix(biDirMinPairs);
tarCong=components(graph);
tarCong=tarCong[[1]];
return(tarCong);
}

### convert group represetations to target set list representations
TargetSetList=function(tarCong)  ##tarCong is a vector (length=|# of objects|
{
uq=unique(tarCong);
len_uq=length(uq);
targetSetList=as.list(numeric(len_uq)); dim(targetSetList)=c(1,len_uq);
  for(i in 1:len_uq)
  {
  i_posInd=which(tarCong==i);
  targetSetList[[i]]=i_posInd;  
  }
return(targetSetList);
}

### recover targetSetSeq into its equivalence classes
TargetSetSeqRec=function(targetSetSeq)
{
len=length(targetSetSeq);
targetSetSeqRec=as.list(numeric(len));dim(targetSetSeqRec)=c(len,1);
catMatrix=matrix(nrow=2,ncol=length(targetSetSeq[[1]]));
catMatrix[1,]=catMatrix[2,]=1:length(targetSetSeq[[1]]);
targetSetSeqRec[[1]]=catMatrix;
   for(i in 2:(len-1))
   {
   list=targetSetSeq[[i]];
   cc=unlist(list);
   dd=unlist(lapply(targetSetSeq[[i]],length));
   cumdd=c(0,cumsum(dd));
   pos_first=cc[(cumdd[1]+1):cumdd[2]];
   location=match(catMatrix[2,],pos_first);
   value=!is.na(location);
   ind=which(value==1);
   catMatrix[2,ind]=1;
      for(k in 2:(length(cumdd)-1))
      {     
      pos=cc[(cumdd[k]+1):cumdd[k+1]];
      location=match(catMatrix[2,],pos);
      value=!is.na(location);
      ind=which(value==1);
      catMatrix[2,ind]=k;
      }
   pos_last=cc[(cumdd[length(dd)]+1):cumdd[length(dd)+1]];
   location=match(catMatrix[2,],pos_last);
   value=!is.na(location);
   ind=which(value==1);
   catMatrix[2,ind]=length(dd);
   targetSetSeqRec[[i]]=catMatrix;
   }
catMatrix[1,]=1:length(targetSetSeq[[1]]);catMatrix[2,]=1;
targetSetSeqRec[[len]]=catMatrix;
return(targetSetSeqRec);
}


## calculate probability matrix  based on rough sets
ProcessPrbMatrix=function(targetSetSeqRec)
{
leng=length(targetSetSeqRec);
##targetSetSeqRecProb=as.list(numeric(leng));dim(targetSetSeqRecProb)=c(leng,1);
cMatrix=targetSetSeqRec[[1]];
processPrbMatrix=matrix(nrow=leng-1,ncol=dim(cMatrix)[2]);
   for(i in 1:(leng-1))
   {
   catMatrix1=targetSetSeqRec[[i]];
   catMatrix2=targetSetSeqRec[[i+1]];
   row1=catMatrix1[2,];
   row2=catMatrix2[2,];
       for(j in 1:dim(cMatrix)[2])
       {
       numerator=sum(row1==row1[j]);
       denominator=sum(row2==row2[j]);
       prb=numerator/denominator;
       processPrbMatrix[i,j]=prb;
       }
   }
return(processPrbMatrix);
}


## calculate probability distribution
ProDis=function(targetSetSeqRec)
{
leng=length(targetSetSeqRec);
nOBJ=dim(targetSetSeqRec[[1]])[2];
proDis=matrix(nrow=leng-1,ncol=nOBJ);
pr=ProcessPrbMatrix(targetSetSeqRec);
proDis[1,]=pr[1,]/length(targetSetSeq[[2]]);
  for(i in 2:(leng-1))
  {
  num_classes=length(targetSetSeq[[i+1]]);
  proDis[i,]=apply(pr[1:i,],2,prod)/num_classes;
  }
return(proDis);
}


## calculate average entropy for the probability distribution
Entropy=function(proDis)
{
logPro=log(10,proDis);
entropy_elementwise=-proDis*logPro;
entropy=sum(entropy_elementwise);
return(entropy);
}

######################################################################################
###############################  Data Analysis  ######################################  
######################################################################################
## import and prepare data (fixed data)
dataMatrix=read.csv(file.choose());
dataMatrix=do.call(cbind,dataMatrix);


## sampled data
size_samples=100;
sample_domain=10;
nOBJ=28;
num_features=14;
samples_entropy=rep(0,size_samples);
for(m in 1:size_samples)
{
dataMatrix=matrix(sample(1:sample_domain,nOBJ*num_features,replace=TRUE),nOBJ,num_features);
## initial target set list
tarCong=1:nOBJ;
targetSetList=TargetSetList(tarCong);

## initial rough distance matrix
disMatrix=matrix(0,nrow=nOBJ,ncol=nOBJ);
   for(i in 1:nOBJ)
   {
   targetSeti=targetSetList[[i]];
   feaRepi=FeaRep(targetSeti,dataMatrix);
       for(j in 1:nOBJ)
       {
       targetSetj=targetSetList[[j]];
       feaRepj=FeaRep(targetSetj,dataMatrix);
       disij_FeaReps=Dis_FeaReps(feaRepi,feaRepj);
       disMatrix[i,j]=disij_FeaReps;
       }
   }

i=1;targetSetSeq=list();
targetSetSeq[[1]]=targetSetList;

num_row=dim(disMatrix)[1];
while(num_row>1)
{i=i+1;
minPosInd=MinPosInd(disMatrix);
biDirMinPairs=BiDirMinPairs(minPosInd);
tarCong=TarCong(biDirMinPairs);
targetSetList=TargetSetList(tarCong);
targetSetSeq[[i]]=targetSetList;
   disMatrix=matrix(0,nrow=length(targetSetList),ncol=length(targetSetList));
   for(j in 1:length(targetSetList))
   {
   targetSetj=targetSetList[[j]];
   feaRepj=FeaRep(targetSetj,dataMatrix);
       for(k in 1:length(targetSetList))
       {
       targetSetk=targetSetList[[k]];
       feaRepk=FeaRep(targetSetk,dataMatrix);
       disjk_FeaReps=Dis_FeaReps(feaRepj,feaRepk);
       disMatrix[j,k]=disjk_FeaReps;
       }
   }
num_row=dim(disMatrix)[1];
}

targetSetSeqRec=TargetSetSeqRec(targetSetSeq);
proDis=ProDis(targetSetSeqRec);
entropy=Entropy(proDis);
samples_entropy[m]=entropy;
}

p11=samples_entropy;

p10=samples_entropy;

p9=samples_entropy;

p8=samples_entropy;

p7=samples_entropy;

p6=samples_entropy;

p5=samples_entropy;

p4=samples_entropy;

p3=samples_entropy;

p2=samples_entropy;

p1=samples_entropy;




par(mfrow=c(3,2));
plot(p2,main="Entropy under 100 samples ", xlab="Samples (nObj=28, nFea=14)",ylab="Entropy")
plot(p1,main="Entropy under 300 samples ", xlab="Samples (nObj=28, nFea=14)",ylab="Entropy")
plot(p3,main="Entropy under 100 samples ", xlab="Samples (nObj=28, nFea=6)",ylab="Entropy")
plot(p4,main="Entropy under 300 samples ", xlab="Samples (nObj=28, nFea=6)",ylab="Entropy")
plot(p5,main="Entropy under 100 samples ", xlab="Samples (nObj=28, nFea=40)",ylab="Entropy")
plot(p6,main="Entropy under 300 samples ", xlab="Samples (nObj=28, nFea=40)",ylab="Entropy")




par(mfrow=c(3,2));
plot(p7,main="Entropy under 100 samples ", xlab="Samples (nObj=15, nFea=14)",ylab="Entropy")
plot(p8,main="Entropy under 300 samples ", xlab="Samples (nObj=15, nFea=14)",ylab="Entropy")
plot(p9,main="Entropy under 100 samples ", xlab="Samples (nObj=60, nFea=14)",ylab="Entropy")
plot(p10,main="Entropy under 300 samples ", xlab="Samples (nObj=60, nFea=14)",ylab="Entropy")









### directional minimal pairs
minPosInd=MinPosInd(disMatrix);

## bi-diretional minimal pairs
biDirMinPairs=BiDirMinPairs(minPosInd);

## targets being conjugated
tarCong=TarCong(biDirMinPairs);
targetSetList=TargetSetList(tarCong);

## rough distance matrix
disMatrix=matrix(0,nrow=length(targetSetList),ncol=length(targetSetList));
for(i in 1:length(targetSetList))
{
targetSeti=targetSetList[[i]];
feaRepi=FeaRep(targetSeti,dataMatrix);
  for(j in 1:length(targetSetList))
  {
  targetSetj=targetSetList[[j]];
  feaRepj=FeaRep(targetSetj,dataMatrix);
  disij_FeaReps=Dis_FeaReps(feaRepi,feaRepj);
  disMatrix[i,j]=disij_FeaReps;
  }
}























##plotting
g=graph_from_adjacency_matrix(minPosInd);

h=graph_from_adjacency_matrix(biDirMinPairs);

plot(g)
plot(h)
graph=graph_from_adjacency_matrix(biDirMinPairs);
integration=components(graph);

feaRep=FeaRep(c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),dataMatrix) 

##testing

###testing
dataColVec=dataMatrix[,1];
eqClass=EqClass(dataColVec);
targetSet=c(1,3,4,5);
lowUpperVecs=LowUpperVecs(targetSet,eqClass);
lowUpperVecs;
bp=BiDirMinPairs(mp)
TarCong(mp)

targetSet=c(1,4);
feaRep=FeaRep(targetSet,dataMatrix);




Smp(disMat)

feaRep1=FeaRep(targetColVec1,dataMatrix);
feaRep2=FeaRep(targetColVec2,dataMatrix);
feaRep3=FeaRep(targetColVec3,dataMatrix);

dis12_FeaReps=Dis_FeaReps(feaRep1,feaRep2);
dis12_FeaReps
dis13_FeaReps=Dis_FeaReps(feaRep1,feaRep3);
dis13_FeaReps
dis23_FeaReps=Dis_FeaReps(feaRep2,feaRep3);
dis23_FeaReps




dataColVec

eqClass=EqClass(dataColVec)

eqClass

lowUpperVecs=LowUpperVecs(targetColVec,eqClass) 

lowUpperVecs


proDis=ProDis(targetSetSeqRec)
entropy=Entropy(proDis);




v1=pr[1,];
v2=apply(pr[1:2,],2,prod);
v3=apply(pr[1:3,],2,prod);
v4=apply(pr[1:4,],2,prod);
v5=apply(pr[1:5,],2,prod);
v6=apply(pr[1:6,],2,prod);
v7=apply(pr[1:7,],2,prod);
v8=apply(pr[1:8,],2,prod);
v9=apply(pr[1:9,],2,prod);


sum(v1);
sum(v2);
sum(v3);
sum(v4);
sum(v5);
sum(v6);
sum(v7);
sum(v8);
sum(v9);


