setwd("/home/invitado/Escritorio/Ale/capitulo2/EN")

DD=read.csv("DDsinBulten.csv")
BBs1=read.csv("BBs1sinBulten.csv")
BBs2=read.csv("BBs2sinBulten.csv")
BBs3=read.csv("BBs3sinBulten.csv")
alfas1=read.csv("alfas1sinBulten.csv")

BB1=as.matrix(BBs1[,2:ncol(BBs1)])
rownames(BB1)=BBs1[,1]
BB2=as.matrix(BBs2[,2:ncol(BBs2)])
rownames(BB2)=BBs2[,1]
BB3=as.matrix(BBs3[,2:ncol(BBs3)])
rownames(BB3)=BBs3[,1]

BB1=BB1[,8:14]
BB2=BB2[,8:14]
BB3=BB3[,8:14]
BB2=1/(1+exp(-BB2))-.5
BB3=-1/(1+exp(-BB3))*5/1000

alphas1=as.matrix(alfas1[,2:ncol(alfas1)])
alphasi=as.matrix(diag(alphas1))
rownames(alphasi)=alfas1[,1]

spnum=nrow(BB1)

fundamentalran=function(DD,BB1,BB2,BB3,alphasi,depth){
	resul=alphasi
	BB=exp(BB1+BB2*depth+BB3*depth^2)
	year=floor(runif(1)*7)+1
	resul=(BB[,year]-DD[,2])/(alphasi*DD[,2])
	return(resul)
}

aaa=fundamentalran(DD,BB1,BB2,BB3,alphasi,10)


#pero como esto no es una simulación usamos la media geométrica de BB??
fundamentalgm=function(DD,BB1,BB2,BB3,alphasi,depth){
	resul=alphasi
	BB=exp(BB1+BB2*depth+BB3*depth^2)
	BB=rowMeans(BB)
	resul=(BB-DD[,2])/(alphasi*DD[,2])
	return(resul)
}

bbb=fundamentalgm(DD,BB1,BB2,BB3,alphasi,5)

fundamentalprof=function(DD,BB1,BB2,BB3,alphasi){
	prof=seq(3,28,0.1) 
	resul=matrix(nrow=nrow(alphasi),ncol=length(prof))
	for(i in 1:length(prof)) resul[,i]=fundamentalgm(DD,BB1,BB2,BB3,alphasi,prof[i])
	return(resul)
}

prue=fundamentalprof(DD,BB1,BB2,BB3,alphasi)
saveRDS(prue,file="Pruesininter")

prue1=prue[-c(5,11,12,14,29,31,32,34,35),]

#Plot Vergnon Function
prof=seq(3,28,0.1) 
#pdf(file="FundamentalSpAbu")
plot(prof,colSums(prue1),type="l",ylab="Soil Depth (cm)",xlab="Species Abundance")
#dev.off()


plot(prof,prue[1,],type="l",ylim=c(0,max(prue)))
for(i in 2:36) lines(prof,prue[i,],col=i)


segder=function(objver){
	sal=2:251
	for(i in 2:251) sal[i-1]=((objver[i+1]-objver[i])/0.004-(objver[i]-objver[i-1])/0.004)/0.004
	sal
}

#Nota: No hay que estandarizar a 1 bajo la curva
StAbu=prue[36,]
StProf=seq(3,28,0.1)

#plot(StProf,StAbu,xlab="Standardized Soil Depth",ylab="Standardized Abundance")
deriv1=segder(StAbu)

saveRDS(StAbu,file="Lam36")
saveRDS(deriv1,file="Lam36deriv")


#read the GAM model for each species
mod1=readRDS("Lam1")
deriv1=readRDS("Lam1deriv")
mod2=readRDS("Lam2")
deriv2=readRDS("Lam2deriv")
mod3=readRDS("Lam3")
deriv3=readRDS("Lam3deriv")
mod4=readRDS("Lam4")
deriv4=readRDS("Lam4deriv")
mod6=readRDS("Lam6")
deriv6=readRDS("Lam6deriv")
mod7=readRDS("Lam7")
deriv7=readRDS("Lam7deriv")
mod8=readRDS("Lam8")
deriv8=readRDS("Lam8deriv")
mod9=readRDS("Lam9")
deriv9=readRDS("Lam9deriv")
mod10=readRDS("Lam10")
deriv10=readRDS("Lam10deriv")
mod13=readRDS("Lam13")
deriv13=readRDS("Lam13deriv")
mod15=readRDS("Lam15")
deriv15=readRDS("Lam15deriv")
mod16=readRDS("Lam16")
deriv16=readRDS("Lam16deriv")
mod17=readRDS("Lam17")
deriv17=readRDS("Lam17deriv")
mod18=readRDS("Lam18")
deriv18=readRDS("Lam18deriv")
mod19=readRDS("Lam19")
deriv19=readRDS("Lam19deriv")
mod20=readRDS("Lam20")
deriv20=readRDS("Lam20deriv")
mod21=readRDS("Lam21")
deriv21=readRDS("Lam21deriv")
mod22=readRDS("Lam22")
deriv22=readRDS("Lam22deriv")
mod23=readRDS("Lam23")
deriv23=readRDS("Lam23deriv")
mod24=readRDS("Lam24")
deriv24=readRDS("Lam24deriv")
mod25=readRDS("Lam25")
deriv25=readRDS("Lam25deriv")
mod26=readRDS("Lam26")
deriv26=readRDS("Lam26deriv")
mod27=readRDS("Lam27")
deriv27=readRDS("Lam27deriv")
mod28=readRDS("Lam28")
deriv28=readRDS("Lam28deriv")
mod30=readRDS("Lam30")
deriv30=readRDS("Lam30deriv")
mod33=readRDS("Lam33")
deriv33=readRDS("Lam33deriv")
mod36=readRDS("Lam36")
deriv36=readRDS("Lam36deriv")

prof=seq(3,28,0.1)
#scale soil depth between 0 and 1

#todas sp
ds=list(deriv1,deriv2,deriv3,deriv4,deriv6,deriv7,deriv8,deriv9,deriv10,deriv13,deriv15,deriv16,deriv17,deriv18,deriv19,deriv20,deriv21,deriv22,deriv23,deriv24,deriv25,deriv26,deriv27,deriv28,deriv30,deriv33,deriv36)
#sin planiv rictri
#ds=list(deriv1,deriv2,deriv3,deriv4,deriv8,deriv9,deriv10,deriv13,deriv14,deriv15,deriv16,deriv17,deriv18,deriv19,deriv20,deriv22,deriv24,deriv25,deriv26,deriv27,deriv29,deriv30,deriv34,deriv36)
#sin tagmic hilcen
#ds=list(deriv1,deriv2,deriv3,deriv4,deriv8,deriv9,deriv10,deriv13,deriv14,deriv15,deriv16,deriv17,deriv18,deriv19,deriv20,deriv21,deriv22,deriv23,deriv24,deriv25,deriv26,deriv27,deriv30,deriv34)

ys=list(mod1,mod2,mod3,mod4,mod6,mod7,mod8,mod9,mod10,mod13,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod30,mod33,mod36)

#Para ver qué tan equitativas están las especies
###########################
todas=c(deriv1,deriv2,deriv3,deriv4,deriv8,deriv9,deriv10,deriv13,deriv15,deriv16,deriv17,deriv18,deriv19,deriv20,deriv21,deriv22,deriv23,deriv24,deriv25,deriv26,deriv27,deriv29,deriv30,deriv34,deriv36)
todasdf=data.frame(sort(rep(1:26,250)),todas)
names(todasdf)=c("sp","segd")
boxplot(segd~sp,data=todasdf)

todas2=data.frame(deriv1,deriv2,deriv3,deriv4,deriv8,deriv9,deriv10,deriv13,deriv14,deriv15,deriv16,deriv17,deriv18,deriv19,deriv20,deriv21,deriv22,deriv23,deriv24,deriv25,deriv26,deriv27,deriv29,deriv30,deriv34,deriv36)
todas2=todas2[-250,]
sdsd=apply(todas2,2,sd)
suma=sdsd/sum(sdsd)
shannon=-sum(suma*log(suma))
shannon

sumabs=function(x) sum(abs(x))
sdsd=apply(todas2,2,sumabs)
suma=sdsd/sum(sdsd)
shannon=-sum(suma*log(suma))
shannon
##############################

signo=function(x) sign(x)*sqrt(abs(x))


derivadas=function(ds){
	res=ds
	for(i in 1:length(ds)){
		ds[[i]]=ds[[i]][-250]
		res[[i]]=signo(ds[[i]])
	}
	res
}

tderivs=derivadas(ds)
#tderivs=derivadas(ds2deriv)
#summary(segder(mod4SSt))
summary(tderivs[[4]])
summary(segder(mod4))



#todas sp
sumd=tderivs[[1]]+tderivs[[2]]+tderivs[[3]]+tderivs[[4]]+tderivs[[5]]+tderivs[[6]]+tderivs[[7]]+tderivs[[8]]+tderivs[[9]]+tderivs[[10]]+tderivs[[11]]+tderivs[[12]]+tderivs[[13]]+tderivs[[14]]+tderivs[[15]]+tderivs[[16]]+tderivs[[17]]+tderivs[[18]]+tderivs[[19]]+tderivs[[20]]+tderivs[[21]]+tderivs[[22]]+tderivs[[23]]+tderivs[[24]]+tderivs[[25]]+tderivs[[26]]+tderivs[[27]]
#sin una sp
#sumd=tderivs[[1]]+tderivs[[2]]+tderivs[[3]]+tderivs[[4]]+tderivs[[5]]+tderivs[[6]]+tderivs[[7]]+tderivs[[8]]+tderivs[[9]]+tderivs[[10]]+tderivs[[11]]+tderivs[[12]]+tderivs[[13]]+tderivs[[14]]+tderivs[[15]]+tderivs[[16]]+tderivs[[17]]+tderivs[[18]]+tderivs[[19]]+tderivs[[20]]+tderivs[[21]]+tderivs[[22]]+tderivs[[24]]+tderivs[[25]]+tderivs[[26]]
#sin dos sp
#sumd=tderivs[[1]]+tderivs[[2]]+tderivs[[3]]+tderivs[[4]]+tderivs[[5]]+tderivs[[6]]+tderivs[[7]]+tderivs[[8]]+tderivs[[9]]+tderivs[[10]]+tderivs[[11]]+tderivs[[12]]+tderivs[[13]]+tderivs[[14]]+tderivs[[15]]+tderivs[[16]]+tderivs[[17]]+tderivs[[18]]+tderivs[[19]]+tderivs[[20]]+tderivs[[21]]+tderivs[[22]]+tderivs[[23]]+tderivs[[24]]
plot(prof[c(-251,-240)],sumd,type="l")

cutpasteder=function(dlist,x){ #ys puede ser ylist
	x2=x[-251]
	x2=x2[-250]
	ysum=rep(0,(length(x2)))
	for(i in 1:length(dlist)){
		original=cbind(x2,dlist[[i]])
		cutx=round(runif(1,1,length(x2)-1))
		new=matrix(nrow=nrow(original),ncol=ncol(original))
		new[1:((nrow(original)-cutx)+1),]=original[cutx:nrow(original),]
		if(cutx==1){
			ysum=ysum+new[,2]
		}
		else{
			new[((nrow(original)-cutx)+2):nrow(original),]=original[1:cutx-1,]
			ysum=ysum+new[,2]
		}
	}
	return(ysum)
}

newver=cutpasteder(tderivs,prof)
x2=prof[-251]
x2=x2[-250]
plot(x2,newver,type="l")

null2=function(dlist,x,iter){
	vernul=matrix(ncol=iter,nrow=length(x)-2)
	for(i in 1:iter){
		vernul[,i]=cutpasteder(dlist,x)
	}
	return(vernul)
}

prue=null2(tderivs,prof,100000)	
#saveRDS(prue,file="null2MAsqrtssinplanivrictri")
#saveRDS(prue,file="null2MAsinhilcentagmic")
#conhilcentagmic
#readRDS('/home/invitado/Escritorio/Ale/capitulo2/EN/Sin Boupol/con inter/100000null2')->prue
#sinhilcentagmic
#readRDS('/home/invitado/Escritorio/Ale/capitulo2/EN/Sin Boupol/con inter/null2sinhilcentagmic')->prue

plotintervals3=function(nullver){
	resul=matrix(nrow=nrow(nullver),ncol=ncol(nullver))
	for(i in 1:nrow(nullver)){
		resul[i,]=sort(nullver[i,])
	}
	return(resul)
}

ord=plotintervals3(prue)

#pdf(file="null2MAsqrtssinplanivrictri.pdf")
pdf(file="SecDerivFundamental.pdf")
plot(prof[c(-251,-240)],ord[,2500],type="l",ylim=c(-60,-10),ylab="Second Derivative",xlab="Soil Depth (cm)")
lines(prof[c(-251,-240)],ord[,97500])
#lines(prof[c(-251,-240)],ord[,95000],col="green")
#lines(prof[c(-251,-240)],ord[,5000],col="green")

lines(prof[c(-251,-240)],sumd,col="red")
dev.off()

#Hellinger Distances
#data should be standardized 
prue1[which(prue1<0)]=0
dat=matrix(nrow=nrow(prue1),ncol=ncol(prue1))
for (i in 1:nrow(prue1)){
	Abusum=sum(prue1[i,])
	dat[i,]=prue1[i,]/Abusum
	} 

hellmat=function(dat){
	res=matrix(nrow=nrow(dat),ncol=nrow(dat))
	for(i in 1:nrow(dat)){
		for(j in 1:nrow(dat)){
			res[i,j]=sqrt(1-sum(sqrt(dat[i,]*dat[j,])))
		}
	}
	res
}

res2=hellmat(dat)
#write.csv(res2,file="hellingerfun.csv")
read.csv("hellingerfunnom.csv")->res2
res=res2[,-1]
rownames(res)=res2[,1]
colnames(res)=res2[,1]

hell=as.dist(res)
hclust(hell,method="ward.D")
plot(hclust(hell,method="ward.D"))
 

#para sacar medias de las especies
dat=read.csv('/home/invitado/Escritorio/Ale/capitulo2/EN/Sin Boupol/datosfunnom.csv')
dat=dat[,-1]

pm=function(x) sum(prof*x)

apply(dat,2,pm)->abumeans
abumeans=sort(abumeans)