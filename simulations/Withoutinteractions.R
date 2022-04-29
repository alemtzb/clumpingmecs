setwd('/home/invitado/Escritorio/Ale/github repo/clumpingmecs/simulations' )

#read data
DD=read.csv("DDs.csv") #mortality parameter
BBs1=read.csv("BBs1.csv") #parameter a
BBs2=read.csv("BBs2.csv") #parameter b
BBs3=read.csv("BBs3.csv") #parametro c
alfas1=read.csv("alfas1.csv") #parameter g 

#convert all data to matrices
BB1=as.matrix(BBs1[,2:ncol(BBs1)])
rownames(BB1)=BBs1[,1]
BB2=as.matrix(BBs2[,2:ncol(BBs2)])
rownames(BB2)=BBs2[,1]
BB3=as.matrix(BBs3[,2:ncol(BBs3)])
rownames(BB3)=BBs3[,1]
alphas1=as.matrix(alfas1[,2:ncol(alfas1)])

#extract intraspecific interactions from the parameter g matrix
alphasi=as.matrix(diag(alphas1))
rownames(alphasi)=alfas1[,1]

#We take into account only the birth rate parameters of the last seven years
BB1=BB1[,8:14]
BB2=BB2[,8:14]
BB3=BB3[,8:14]

#Detransform parameters
BB2=1/(1+exp(-BB2))-.5
BB3=-1/(1+exp(-BB3))*5/1000

#To obtain total abundance of the species  for which we did not calculate pairwise interactions
spnum=ncol(alphas1) #to obtain the number of species in our study

#Function to simulate the abundance of species in the next year
fundamentalgm=function(DD,BB1,BB2,BB3,alphasi,depth){
	resul=alphasi
	BB=exp(BB1+BB2*depth+BB3*depth^2)
	BB=rowMeans(BB)
	resul=(BB-DD[,2])/(alphasi*DD[,2])
	return(resul)
}

#Function that runs the simulation 
fundamentalprof=function(DD,BB1,BB2,BB3,alphasi){
	prof=seq(3,28,0.1) 
	resul=matrix(nrow=nrow(alphasi),ncol=length(prof))
	for(i in 1:length(prof)) resul[,i]=fundamentalgm(DD,BB1,BB2,BB3,alphasi,prof[i])
	return(resul)
}

prue=fundamentalprof(DD,BB1,BB2,BB3,alphasi)

#To plot the outcome. Eah color represents a different species
prof=seq(3,28,0.1) 
plot(prof,prue[1,],type="l",ylim=c(0,max(prue)))
for(i in 2:36) lines(prof,prue[i,],col=i)


