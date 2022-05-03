#Alejandra Martínez Blancas & Carlos Martorell 03/05/22 alemtzb@ciencias.unam.mx
#Code for simulating species abundances when the other species in their clusters have been removed

#read data
DD=read.csv("clumpingmecs/simulations/DDs.csv") #mortality parameter
BBs1=read.csv("clumpingmecs/simulations/BBs1.csv") #parameter a
BBs2=read.csv("clumpingmecs/simulations/BBs2.csv") #parameter b
BBs3=read.csv("clumpingmecs/simulations/BBs3.csv") #parametro c
alfas1=read.csv("clumpingmecs/simulations/alfas1.csv") #parameter g 
alfas2=read.csv("clumpingmecs/simulations/alfas2.csv") #parameter h
betas=read.csv("clumpingmecs/simulations/betas.csv")#facilitation parameter

#convert all data to matrices
BB1=as.matrix(BBs1[,2:ncol(BBs1)])
rownames(BB1)=BBs1[,1]
BB2=as.matrix(BBs2[,2:ncol(BBs2)])
rownames(BB2)=BBs2[,1]
BB3=as.matrix(BBs3[,2:ncol(BBs3)])
rownames(BB3)=BBs3[,1]
alphas1=as.matrix(alfas1[,2:ncol(alfas1)])
rownames(alphas1)=alfas1[,1]
alphas2=as.matrix(alfas2[,2:ncol(alfas2)])
rownames(alphas2)=alfas2[,1]
bbetas=as.matrix(betas[,2:ncol(betas)])
rownames(bbetas)=betas[,1]

#We take into account only the birth rate parameters of the last seven years
BB1=BB1[,8:14]
BB2=BB2[,8:14]
BB3=BB3[,8:14]
#Detransform parameters
BB2=1/(1+exp(-BB2))-.5
BB3=-1/(1+exp(-BB3))*5/1000
-999->alphas1[which(is.na(alphas1[,])=="TRUE")]
0->alphas2[which(is.na(alphas2[,])=="TRUE")]
bbetas=1/(1+exp(-bbetas))
0->bbetas[which(is.na(bbetas[,])=="TRUE")]

#To obtain total abundance of the species  for which we did not calculate pairwise interactions
spnum=ncol(alphas1)
tx0=matrix(ncol=ncol(alphas1),nrow=nrow(alphas1))
0->tx0[which(is.na(alphas1[,])=="FALSE")]
1->tx0[which(is.na(alphas1[,])=="TRUE")]
tx0=tx0[,-37]

#separate the parameters of the species for which we have abundance data
alpabu1=alphas1[1:33,1:37]
alpabu2=alphas2[1:33,1:37]
betabu=bbetas[1:33,1:37]
DDabu=DD[1:33,]
BB1abu=BB1[1:33,]
BB2abu=BB2[1:33,]
BB3abu=BB3[1:33,]
#separate the parameters of the species for which we have presence/absence data
alppa1=alphas1[34:36,1:37]
alppa2=alphas2[34:36,1:37]
betpa=bbetas[34:36,1:37]
DDpa=DD[34:36,]
BB1pa=BB1[34:36,]
BB2pa=BB2[34:36,]
BB3pa=BB3[34:36,]

#Function for species with abundance data to simulate the abundance of species in the next year
lam=function(DD,BB1,BB2,BB3,alphas1,alphas2,bbetas,tx,txothers,year,depth){
	surv=(1-DD[,2])*tx[1:33]
	BB1=BB1[,year]
	BB2=BB2[,year]
	BB3=BB3[,year]
	BB=exp(BB1+BB2*depth+BB3*depth^2)
	alphas=exp(alphas1+alphas2*depth)
	alphaspre=alphas[,1:36]%*%tx 
	alphasothers=alphas[,37]*txothers
	alphasum=alphaspre+alphasothers
	txmas=log(tx+1)
	txmasother=log(txothers+1)
	betaspre=bbetas[,1:36]%*%txmas
	betasothers=bbetas[,37]*txmasother
	betassum=betaspre+betasothers
	fac=exp(betassum)
	new=BB*tx[1:33]/(1+alphasum)*fac
	t2=surv+new
	return(t2)
}


#Function for species with presence absence daata to simulate if the species is present of absent in the next year
lampa=function(DD,BB1,BB2,BB3,alphas1,alphas2,bbetas,tx,txothers,year,depth){
	surv=(1-DD[,2])*tx[34:36]
	BB1=BB1[,year]
	BB2=BB2[,year]
	BB3=BB3[,year]
	BB=exp(BB1+BB2*depth+BB3*depth^2)
	alphas=exp(alphas1+alphas2*depth) 
	alphaspre=alphas[,1:36]%*%tx 
	alphasothers=alphas[,37]*txothers
	alphasum=alphaspre+alphasothers
	txmas=log(tx+1)
	txmasother=log(txothers+1)
	betaspre=bbetas[,1:36]%*%txmas
	betasothers=bbetas[,37]*txmasother
	betassum=betaspre+betasothers
	fac=exp(betassum)
	new=BB*tx[34:36]/(1+alphasum)*fac
	t2=1-(1-surv)*exp(-new)
	return(t2)
}



#To obtain initial abuncances of the species
spnum=ncol(alphas1)
tx=matrix(ncol=1,nrow=nrow(alphas1))
tx[1:33,1]=runif(33, min = .001, max = 1)
tx[34:36,1]=runif(3, min = .001, max = 1)
tx[35,1]=0 #eliminates species 35 which we weren't able to model correctly
rownames(tx)=rownames(alphas1)


#Create data where each species is present without the other species in their cluster
gps=matrix(tx,nrow=nrow(tx),ncol=36)
rownames(gps)=rownames(tx)
gps[c(20,22,29,35),13]=0 #Hetpin stays
gps[c(13,22,29,35),20]=0 #Pheoli stays
gps[c(13,20,29,35),22]=0 #Porpil stays
gps[c(13,20,22,35),29]=0 #Tagmic stays
gps[c(2,3,4,8,10,14,15,16,17,19,23,24,25,26,30,35),1]=0 #Ariads stays
gps[c(1,3,4,8,10,14,15,16,17,19,23,24,25,26,30,35),2]=0 #Aridiv
gps[c(1,2,4,8,10,14,15,16,17,19,23,24,25,26,30,35),3]=0 #Bouhir
gps[c(1,2,3,8,10,14,15,16,17,19,23,24,25,26,30,35),4]=0 #Bousco
gps[c(1,2,3,4,10,14,15,16,17,19,23,24,25,26,30,35),8]=0 #Eupmen
gps[c(1,2,3,4,8,14,15,16,17,19,23,24,25,26,30,35),10]=0 #Fĺoped
gps[c(1,2,3,4,8,10,15,16,17,19,23,24,25,26,30,35),14]=0 #Mickun
gps[c(1,2,3,4,8,10,14,16,17,19,23,24,25,26,30,35),15]=0 #Milbif
gps[c(1,2,3,4,8,10,14,15,17,19,23,24,25,26,30,35),16]=0 #Muhper
gps[c(1,2,3,4,8,10,14,15,16,19,23,24,25,26,30,35),17]=0 #Muhpha
gps[c(1,2,3,4,8,10,14,15,16,17,23,24,25,26,30,35),19]=0 #Oxalun
gps[c(1,2,3,4,8,10,14,15,16,17,19,24,25,26,30,35),23]=0 #Rictri
gps[c(1,2,3,4,8,10,14,15,16,17,19,23,25,26,30,35),24]=0 #Sanpro
gps[c(1,2,3,4,8,10,14,15,16,17,19,23,24,26,30,35),25]=0 #Sedote
gps[c(1,2,3,4,8,10,14,15,16,17,19,23,24,25,30,35),26]=0 #Spoten
gps[c(1,2,3,4,8,10,14,15,16,17,19,23,24,25,26,35),30]=0 #Thyaur
gps[c(18,21,27,34,35),9]=0 #Evoser
gps[c(9,21,27,34,35),18]=0 #Muhrig
gps[c(9,18,27,34,35),21]=0 #Planiv
gps[c(9,18,21,34,35),27]=0 #Stedul
gps[c(9,18,21,27,35),34]=0 #Boucho
gps[c(36,35),32]=0 #Tripur
gps[c(32,35),36]=0 #Hilcen

#Function that integrates both species with abuncance data and species with presence absence data into the same simulation
lamx=function(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,year,depth){
	txothers=as.vector(tx0%*%tx)
	t2abu=lam(DDabu,BB1abu,BB2abu,BB3abu,alpabu1,alpabu2,betabu,tx,txothers[1:33],year,depth)
	t2pa=lampa(DDpa,BB1pa,BB2pa,BB3pa,alppa1,alppa2,betpa,tx,txothers[34:36],year,depth)
	t2=matrix(nrow=36,ncol=1)
	t2[1:33,]=t2abu[1:33,]
	t2[34:36,]=t2pa
	rownames(t2)=rownames(DD)
	return(t2)
}

#Function that runs the simulation for a species number of time steps (iter) but ilimanated a specific number of runs (burn)
simu=function(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,depth,iter=1000,burn=100){
	for(i in 1:burn) {
		tx=lamx(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,floor(runif(1)*7)+1,depth)
	}
	sal=matrix(nrow=36,ncol=iter)
	for(i in 1:iter) {
		sal[,i]=tx
		tx=lamx(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,floor(runif(1)*7)+1,depth)
	}


		rownames(sal)=rownames(BB1)
		return(sal)
}

#Runs the simulation for each soil depth between 3 and 28 cm every 0.1 cm
profs=function(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,iter){
	prof=seq(3,28,0.5)
	sal=matrix(nrow=36,ncol=length(prof))
	for(i in 1:length(prof)){
		sal[,i]=rowMeans(simu(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,prof[i],iter=iter))
	}
	return(sal)
}

#Runs the simulation where each species is without the other specie in its cluster along the entire soil depth gradient
gscen=function(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx0,iter,gps){
	prof=seq(3,28,0.5)
	resul=array(dim=c(nrow(gps),ncol=length(prof),ncol(gps)))
	rownames(resul)=rownames(alphas1)
	for(i in 1:ncol(gps)) resul[,,i]=profs(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,gps[,i],tx0,iter)
	return(resul)
}


esce=gscen(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx0,100000,gps)

