setwd('/home/invitado/Escritorio/Ale/github repo/clumpingmecs/simulations' )

#read data
DD=read.csv("DDs.csv") #mortality parameter
BBs1=read.csv("BBs1.csv") #parameter a
BBs2=read.csv("BBs2.csv") #parameter b
BBs3=read.csv("BBs3.csv") #parametro c
alfas1=read.csv("alfas1.csv") #parameter g 
alfas2=read.csv("alfas2.csv") #parameter h
betas=read.csv("betas.csv")#facilitation parameter which will later be ignored in the simulation

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
spnum=ncol(alphas1) #to obtain the number of species in our study
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
	new=BB*tx[1:33]/(1+alphasum) #note that facilitation is not included in the calculation 
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
	new=BB*tx[34:36]/(1+alphasum) #note that facilitation is not included in the calculation 
	t2=1-(1-surv)*exp(-new)
	return(t2)
}


#To obtain initial abuncances of the species
spnum=ncol(alphas1)
tx=matrix(ncol=1,nrow=nrow(alphas1))
tx[1:33,1]=runif(33, min = .001, max = 1)
tx[34:36,1]=runif(3, min = .001, max = 1)
tx[35,1]=0 #eliminates species 35

#Function that integrates both species with abuncance data and species with presence absence data into the same simulation
lamx=function(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,year,depth){
	#species with abundance data
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
	prof=seq(3,28,0.1)
	ncat=length(seq(3,28,0.1))
	sal=matrix(nrow=36,ncol=ncat)
	for(i in 1:ncat){
		sal[,i]=rowMeans(simu(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,prof[i],iter=iter))
	}
	return(sal)
}

prue=profs(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,500000)

#To plot the outcome. Eah color represents a different species
plot(-1000,-1000,xlim=c(0,28),ylim=c(0,max(prue)))
prof=seq(3,28,0.1) 
for(i in 1:36) lines(prof,prue[i,],col=i)

#Runs the simulation along the soil depht a specific number of times (rep) and averages the outcome to obtain smoother abuncance curves over the gradient
meanprof=function(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,iter,rep){
	prof=seq(3,28,0.1)
	ncat=length(seq(3,28,0.1))
	sal=array(dim=c(36,ncat,rep))
	for(k in 1:rep){
		sal[,,k]=profs(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,iter)
	}
	return(sal)
}

pruemean=meanprof(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,10000,20)

#To plot the outcome
plot(-1000,-1000,xlim=c(0,28),ylim=c(0,max(pruemean)),xlab="Profundidad de Suelo",ylab="Abundancia")
prof=seq(3,28,0.1) 
for(i in 1:36) lines(prof,rowMeans(pruemean[i,,]),col=i)