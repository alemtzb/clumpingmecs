setwd("/home/invitado/Escritorio/Ale/capitulo2/EN")
#setwd("/home/invitado/Escritorio/Ale/capitulo2/EN/Sin Boupol")

DD=read.csv("DDsinBulten.csv")
BBs1=read.csv("BBs1sinBulten.csv")
BBs2=read.csv("BBs2sinBulten.csv")
BBs3=read.csv("BBs3sinBulten.csv")
alfas1=read.csv("alfas1sinBulten.csv")
alfas2=read.csv("alfas2sinBulten.csv")
betas=read.csv("betassinBulten.csv")
Thetas=read.csv("Thetas.csv")

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

BB1=BB1[,8:14]
BB2=BB2[,8:14]
BB3=BB3[,8:14]
BB2=1/(1+exp(-BB2))-.5
BB3=-1/(1+exp(-BB3))*5/1000
-999->alphas1[which(is.na(alphas1[,])=="TRUE")]
0->alphas2[which(is.na(alphas2[,])=="TRUE")]
bbetas=1/(1+exp(-bbetas))
0->bbetas[which(is.na(bbetas[,])=="TRUE")]

#obtain total abundance of "other" species
spnum=ncol(alphas1)
#tx=matrix(ncol=1,nrow=ncol(alphas1))
#tx[,1]=rep(.01,spnum)
tx0=matrix(ncol=ncol(alphas1),nrow=nrow(alphas1))
0->tx0[which(is.na(alphas1[,])=="FALSE")]
1->tx0[which(is.na(alphas1[,])=="TRUE")]
tx0=tx0[,-37]

#species with abundance data
alpabu1=alphas1[1:33,1:37]
alpabu2=alphas2[1:33,1:37]
betabu=bbetas[1:33,1:37]
DDabu=DD[1:33,]
BB1abu=BB1[1:33,]
BB2abu=BB2[1:33,]
BB3abu=BB3[1:33,]
#species with presence/absence data
alppa1=alphas1[34:36,1:37]
alppa2=alphas2[34:36,1:37]
betpa=bbetas[34:36,1:37]
DDpa=DD[34:36,]
BB1pa=BB1[34:36,]
BB2pa=BB2[34:36,]
BB3pa=BB3[34:36,]




lam=function(DD,BB1,BB2,BB3,alphas1,alphas2,bbetas,tx,txothers,year,depth){
	surv=(1-DD[,2])*tx[1:33]
	BB1=BB1[,year]
	BB2=BB2[,year]
	BB3=BB3[,year]
	BB=exp(BB1+BB2*depth+BB3*depth^2)
	alphas=exp(alphas1+alphas2*depth)
	alphaspre=alphas[,1:36]%*%tx #revisar este paso
	alphasothers=alphas[,37]*txothers
	alphasum=alphaspre+alphasothers
	#revisar este paso
	txmas=log(tx+1)
	txmasother=log(txothers+1)
	betaspre=bbetas[,1:36]%*%txmas
	betasothers=bbetas[,37]*txmasother
	betassum=betaspre+betasothers
	fac=exp(betassum)
	#finalalphas=alphas*tx
	new=BB*tx[1:33]/(1+alphasum)*fac
	t2=surv+new
	return(t2)
}


#species with presence/absence data

lampa=function(DD,BB1,BB2,BB3,alphas1,alphas2,bbetas,tx,txothers,year,depth){
	surv=(1-DD[,2])*tx[34:36]
	BB1=BB1[,year]
	BB2=BB2[,year]
	BB3=BB3[,year]
	BB=exp(BB1+BB2*depth+BB3*depth^2)
	alphas=exp(alphas1+alphas2*depth) 
	alphaspre=alphas[,1:36]%*%tx #revisar este paso
	alphasothers=alphas[,37]*txothers
	alphasum=alphaspre+alphasothers
 #revisar este paso
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



#t2=matrix(nrow=36,ncol=1)
#t2[1:34,]=t2abu[1:35,]
#t2[34:36,]=t2pa


#Correr desde aquí
spnum=ncol(alphas1)
tx=matrix(ncol=1,nrow=nrow(alphas1))
tx[1:33,1]=runif(33, min = .001, max = 1)
tx[34:36,1]=runif(3, min = .001, max = 1)
#tx[1:33,1]=.000001
#tx[34:36,1]=.000001
#tx[c(1:4,9,10,12:14,16:19,21,23,24,26,27,29,30,34)]=.000001
#tx[c(5:8,11,15,20,22,25,28,31:33,35,36)]=1
tx[35,1]=0



#decirle cuales son los datos de abundancia y de pres aus desde afuera
#other data también tiene que ir desde afuera

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

lamx(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,2,20)


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

simu(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,20,iter=10000)->aa

profs=function(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,iter){
	prof=seq(3,28,0.1)
	ncat=length(seq(3,28,0.1))
	sal=matrix(nrow=36,ncol=ncat)
	for(i in 1:ncat){
		sal[,i]=rowMeans(simu(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,prof[i],iter=iter))
	}
	return(sal)
}

profs(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,100000)->prue
saveRDS(prue,file="ProfssinBoupoliter100000")
plot(-1000,-1000,xlim=c(0,28),ylim=c(0,max(prue)))
prof=seq(3,28,0.1) 
for(i in 1:36) lines(prof,prue[i,],col=i)

#meanminmax=matrix(nrow=36,ncol=3)
#for(i in 1:36){
	#meanminmax[i,1]=mean(prue[i,])
	#meanminmax[i,2]=min(prue[i,])
	#meanminmax[i,3]=max(prue[i,])
#}

meanprof=function(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,iter,rep){
	prof=seq(3,28,0.1)
	ncat=length(seq(3,28,0.1))
	sal=array(dim=c(36,ncat,rep))
	for(k in 1:rep){
		sal[,,k]=profs(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,iter)
	}
	return(sal)
}

#pruemean=meanprof(alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,5000,100)
#saveRDS(pruemean,file="Pruemean100sinBulten")
pruemean=meanprof(DDabu,alpabu1,alpabu2,betabu,BB1abu,BB2abu,BB3abu,alppa1,alppa2,betapa,DDpa,BB1pa,BB2pa,BB3pa,tx,tx0,5000,20)
saveRDS(pruemean,file="Pruemean400sinBoupoliter100000")
plot(-1000,-1000,xlim=c(3,28),ylim=c(0,max(pruemean,na.rm=TRUE)),xlab="Profundidad de Suelo",ylab="Abundancia")
prof=seq(3,28,0.1) 
for(i in 1:36) lines(prof,rowMeans(pruemean[i,,]),col=i)


#graficas artículo
azul=c(13,29,20,22)
verde=c(1,2,3,4,8,10,14,15,16,17,19,23,24,25,26,30)
naranja=c(9,18,21,27,35)
rojo=c(32,36)

setwd("/home/invitado/Escritorio/Ale/capitulo2/EN/figuras")
pruemean=readRDS('/home/invitado/Escritorio/Ale/capitulo2/EN/Sin Boupol/Pruemean400sinBoupoliter100000')

pdf(file="Sp abundance all interactions.pdf")
plot(-1000,-1000,xlim=c(3,28),ylim=c(0,max(pruemean,na.rm=TRUE)),xlab="Soil Depth",ylab="Abundance")
prof=seq(3,28,0.1) 
for(i in 1:length(azul)) lines(prof,rowMeans(pruemean[azul[i],,]),col="blue")
for(i in 1:length(verde)) lines(prof,rowMeans(pruemean[verde[i],,]),col="green")
for(i in 1:length(naranja)) lines(prof,rowMeans(pruemean[naranja[i],,]),col="orange")
for(i in 1:length(rojo)) lines(prof,rowMeans(pruemean[rojo[i],,]),col="firebrick1")
dev.off()

#para estandarizar curvas:
rowMeans(pruemean[naranja[i],,])/sum(rowMeans(pruemean[naranja[i],,]))

#sin fac
pruemean=readRDS('/home/invitado/Escritorio/Ale/capitulo2/EN/Sin Boupol/Pruemean400sinBoupolsinfaciter50000')

azul=29
azulcielo=c(20,22,8)
verdelimon=19 
aqua=c(25,15,10) 
verde=c(1,2,3,4,26,30) 
verdecaca=c(16,17,23,24,27)
naranja=c(9,18,21) 

pdf(file="Sp abundance no fac.pdf")
plot(-1000,-1000,xlim=c(3,28),ylim=c(0,max(pruemean,na.rm=TRUE)),xlab="Soil Depth",ylab="Abundance")
prof=seq(3,28,0.1) 
lines(prof,rowMeans(pruemean[azul,,]),col="blue")
for(i in 1:length(azulcielo)) lines(prof,rowMeans(pruemean[azulcielo[i],,]),col="deepskyblue1")
lines(prof,rowMeans(pruemean[verdelimon,,]),col="green4")
for(i in 1:length(aqua)) lines(prof,rowMeans(pruemean[aqua[i],,]),col="springgreen")
for(i in 1:length(verde)) lines(prof,rowMeans(pruemean[verde[i],,]),col="green")
for(i in 1:length(verdecaca)) lines(prof,rowMeans(pruemean[verdecaca[i],,]),col="olivedrab3")
for(i in 1:length(naranja)) lines(prof,rowMeans(pruemean[naranja[i],,]),col="orange")
dev.off()

#sin inter
prue=readRDS('/home/invitado/Escritorio/Ale/capitulo2/EN/Pruesininter')

pdf(file="Sp abundance no inter.pdf")
prof=seq(3,28,0.1) 
plot(prof,prue[1,],type="l",ylim=c(0,max(prue)),xlab="Soil Depth")
for(i in 2:36) lines(prof,prue[i,])
dev.off()