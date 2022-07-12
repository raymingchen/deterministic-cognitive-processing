########## machine learning the optimal (parameters) cognitive process
## Euclidean Distance Matrix (given objects and features)
rm(list=ls());
EucDisMat=function(ExpData)  ## ExpData:NumObj-by-NumFea
{
nobj=dim(ExpData)[1];
nfe=dim(ExpData)[2];
eucdismat=matrix(nrow=nobj,ncol=nobj);
  for(i in 1:nobj)
  {
    for(j in 1:nobj)
    {
    vec_i=ExpData[i,];
    vec_j=ExpData[j,];
    difvec=vec_i-vec_j;
    eucdismat[i,j]=sqrt(difvec%*%difvec);
    }
  }
return(eucdismat);
}

### find the minima and posisition row-wise for matrix, excluding diag
ROWmin=function(MATRIX)   ##MATRIX:a square matrix 
{
diag(MATRIX)=NA;
nmat=dim(MATRIX)[1];  ## row dimension of the matrix
rowmin=matrix(nrow=2,ncol=nmat);
  for(k in 1:nmat)
  {
  mini=min(MATRIX[k,],na.rm=TRUE);
  pos_mini=which(MATRIX[k,]==mini);
  rowmin[1,k]=mini[1];
  rowmin[2,k]=pos_mini[1];
  }
return(rowmin);
}
## grouping via minimal distance
Group=function(Matrix_2_n)  ## Matrix_2_n:2-by-n matrix,1st row:minima, 2nd row:posistions   
{
uq=unique(Matrix_2_n[1,]);
luq=length(uq);
PT=as.list(numeric(luq)); dim(PT)=c(luq,1);
for(k in 1:luq)
{
PT[[k]]=which(Matrix_2_n[1,]==uq[k]); ##PT:partion
}
return(PT);    
}
#### Detect the position of a list whose element length > 1
PosListL1=function(LIST)  ## LIST:n-by-1
{
llist=length(LIST);
PosLs=vector();
j=0;
for(i in 1:llist)
{
   if(length(LIST[[i]])>1)
   {j=j+1;
   PosLs[j]=i}else
   {
    j=j;
   }
}
return(PosLs);
}
###### Convert the positioned matrix into a vector
vect=function(LIST)
{
MAT=do.call(cbind,LIST[PosListL1(LIST)]);
return(as.vector(MAT));
}
###### Hausdorf Matrix
Hdis=function(Ind1,Ind2,DisMat)  ## Ind1,Ind2:index vectors
{
subDismat=as.matrix(DisMat[Ind1,Ind2]);
rowmin=apply(subDismat,1,min,na.rm=TRUE);
colmin=apply(subDismat,2,min,na.rm=TRUE);
rowMax=max(rowmin);
colMax=max(colmin);
MAX=max(rowMax,colMax);
return(MAX);
}
#######################################################
############### DATA FEEDING  ########################
######################################################
## experimental data creation
nobj=28;        ##nobj:number of objects
nfe=14;         ##nfe:number of features
ExpData=matrix(nrow=nobj,ncol=nfe);   ## ExpData:Experimental Data Created
for(i in 1:nobj)
{
ExpData[i,]=sample(seq(1,40,by=1),nfe, replace=TRUE);
}

#### data analysis
Nodes=list();
distMAT=list();
Rowmin=list();
Gp=list(list());
for(i in 1:nobj)
{
Gp[[1]][i]=i;
}

distMAT[[1]]=EucDisMat(ExpData);

i=1;
while(length(Gp[[i]])>=2)
{
i=i+1;
Rowmin[[i]]=ROWmin(distMAT[[i-1]]);
Gp[[i]]=Group(Rowmin[[i]]); 
nNodes_iplus1=length(Gp[[i]]);
distanceMat=matrix(nrow=nNodes_iplus1,ncol=nNodes_iplus1);
   for(j in 1:nNodes_iplus1)
   {
       for(k in 1:nNodes_iplus1)
       {
       distanceMat[j,k]=Hdis(Gp[[i]][[j]],Gp[[i]][[k]],distMAT[[i-1]]);
       }
   }
distMAT[[i]]=distanceMat;
}

## convert Gp to edges
nLevel=length(Gp); ##nLevel:length of levels
NameNodes=list();
for(j in 1:nLevel)
{
nRowGp_j=length(Gp[[j]]);  ##number of rows in Gp_j
Mat=matrix(nrow=2,);
   for(i in 1:nRowGp_j)
   {
   vec=Gp[[j]][[i]];
   matr=matrix(nrow=2,ncol=length(vec));
   matr[1,]=matrix(i,length(vec));
   matr[2,]=vec;
   Mat=cbind(Mat,matr);
   }
NameNodes[[j]]=Mat[,2:(dim(Mat)[2])];
}

### compute Nodes implicitly
Nodes=list();
for(k in 1:nLevel)
{
vlist=list();
nColNameNodes=length(unique(NameNodes[[k]][1,]));
   for(v in 1:nColNameNodes)
   {  
   Mk=NameNodes[[k]];
   vlist[[v]]=Mk[2,Mk[1,]==v];
   }
Nodes[[k]]=vlist;
}

### compute Nodes explicitly
NodesEx=list(list());
NodesEx[[1]]=Nodes[[1]];
for(k in 2:nLevel)
{
MM=list();
    for(i in 1:length(Nodes[[k]]))
    { 
    M=vector();
        for(j in Nodes[[k]][[i]])
        {
        M=append(M,NodesEx[[k-1]][[j]])
        }
    MM[[i]]=M;
    }
NodesEx[[k]]=MM;
}

######## distance for levelled trees
##construction of levelled trees from Gp and Nodes








