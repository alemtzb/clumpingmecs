#Alejandra Martínez Blancas & Ian Xul 11/04/21 alemtzb@ciencias.unam.mx
require(TMB)
require(lhs)
require(deldir)

setwd("/home/invitado/Escritorio/Ale/capitulo2/Poisson")

cleanDt=function(dtSp){
    #Remove these combinations of site/year:
    # Agua blanca F (3) - 2001,2,3,4,5
    dtSp=dtSp[!(dtSp$Localidad=='Agua_Blanca' & dtSp$Exclusion=='Fuera' & dtSp$Cuadro==3 & is.element(dtSp$year,2001:2005)),]
    # Agua mosca F (2,3,6,7,8) - 1:8
    dtSp=dtSp[!(dtSp$Localidad=='Aguamosca' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(2,3,6,7,8)) & is.element(dtSp$year,2001:2008)),]
    # Agua mosca F (8) - 10
    dtSp=dtSp[!(dtSp$Localidad=='Aguamosca' & dtSp$Exclusion=='Fuera' & dtSp$Cuadro==8 & dtSp$year==2010),]
    # Biznaga F (2) - 1,2
    dtSp=dtSp[!(dtSp$Localidad=='Biznaga' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(2)) & is.element(dtSp$year,c(2001,2002))),]
    # Calle de piedra F (1:8) - 1:11
    dtSp=dtSp[!(dtSp$Localidad=='' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,1:8) & is.element(dtSp$year,c(2001:2008,2010,2011))),]
    # Camino real F (2:8) - 1:11
    dtSp=dtSp[!(dtSp$Localidad=='Camino_Real' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,2:8) & is.element(dtSp$year,c(2001:2008,2010,2011))),]
    # Canada moral F (1,3,4,6,7) - 1:11
    dtSp=dtSp[!(dtSp$Localidad=='Canada_Moral' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(1,2,3,6,7)) & is.element(dtSp$year,c(2001:2008,2010,2011))),]
    # El gavilan F (1:7) - 1:8
    dtSp=dtSp[!(dtSp$Localidad=='El_Gavilan' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(1:7)) & is.element(dtSp$year,c(2001:2008))),]
    # El humo F (1:5,7,8) - 1:11
    dtSp=dtSp[!(dtSp$Localidad=='El_Humo' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(1:5,7,8)) & is.element(dtSp$year,c(2001:2008,2010,2011))),]
    # El ovni F (1,5,6,8) - 1:11
    dtSp=dtSp[!(dtSp$Localidad=='El_Ovni' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(1,5,6,8)) & is.element(dtSp$year,c(2001:2008,2010,2011))),]
    # El tambor F (2,4,5,6) - 1:11
    dtSp=dtSp[!(dtSp$Localidad=='El_Tambor' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(2,4,5,6)) & is.element(dtSp$year,c(2001:2008,2010,2011))),]
    # La chilacayota F (4:8) - 1:8
    dtSp=dtSp[!(dtSp$Localidad=='La_Chilacayota' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(4:8)) & is.element(dtSp$year,c(2001:2008))),]
    # Llano de la estrella F (1,2,4,5,6,7) - 1:8
    dtSp=dtSp[!(dtSp$Localidad=='Llano_de_la_estrella' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(1,2,4,5,6,7)) & is.element(dtSp$year,c(2001:2008))),]
    # Loma del coyote F (1:7) - 1:11 
    dtSp=dtSp[!(dtSp$Localidad=='Loma_del_coyote' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(1:7)) & is.element(dtSp$year,c(2001:2008,2010,2011))),]
    # Nachininge F (1,2,4,5,8,10,11) - 1:8
    dtSp=dtSp[!(dtSp$Localidad=='Nachiningue' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(1,2,4,5,8,10,11)) & is.element(dtSp$year,c(2001:2008))),]
    # Nadenda F (2,4:8) - 1:8
    dtSp=dtSp[!(dtSp$Localidad=='Nadenda' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(2,4:8)) & is.element(dtSp$year,c(2001:2008))),]
    # Pedrera cima F (1,3:7) - 1:8
    dtSp=dtSp[!(dtSp$Localidad=='Pedrera_cima' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(1,3:7)) & is.element(dtSp$year,c(2001:2008))),]
    # Pena F (8) - 10,11,12
    dtSp=dtSp[!(dtSp$Localidad=='Pena' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(8)) & is.element(dtSp$year,c(2010,2011,2012))),]
    # Piedradura F (1,3:8) - 1:11
    dtSp=dtSp[!(dtSp$Localidad=='Piedradura' & dtSp$Exclusion=='Fuera' & is.element(dtSp$Cuadro,c(1,3:8)) & is.element(dtSp$year,c(2001:2008,2010,2011))),]

return(dtSp)
}


cleanzero=function(dtSp){
    # Remove entries for squares where the whole square did not have a single individual    
	numcuad=nrow(dtSp)/20
	elim=c()
	for(i in 1:numcuad){
			if(sum(dtSp$t1[((i-1)*20+1):(i*20)])==0) elim=c(elim,((i-1)*20+1):(i*20))
	}
	return(dtSp[-elim,])
}

setupVals = function(dtSp){
    #squares list will hold the different combinations of factors that determine the squares
    squares=list()
    # for loop counts over all unique conbinations of factors, assigning a unique id (count) to each
    count=1
    for(l in unique(dtSp$Localidad)){
        for(e in unique(dtSp$Exclusion[which(dtSp$Localidad==l)])){
            for(c in unique(dtSp$Cuadro[which(dtSp$Localidad==l & dtSp$Exclusion==e)])){
                squares=rbind(squares,list(count,l,e,c))
                count=count+1
            }
        }
    }

    #messSq contains the coordinates of each subsquare in each posible square
    messSq=array(NaN,c(2,20,length(squares[,1])))

    for(i in 1:length(dtSp[,1])){
        #isquare is the index of the square in the squares list which matches the factors of the ith element of dtSp
        isquare=which(squares[,2]==dtSp[i,]$Localidad & squares[,3]==dtSp[i,]$Exclusion & squares[,4]==dtSp[i,]$Cuadro)
        
        # the the data are filled in in order, and they are no longer added when the list is full. There might be faster ways to do this but for now this is fine.
        coords=c(dtSp[i,]$x,dtSp[i,]$y)
        if(sum(is.nan(messSq[,,isquare]))!=0){
            messSq[,,isquare][,which(is.nan(messSq[,,isquare][1,]))[1]]=coords
        }
    }

    # ZZZ holds n 100x20 matrices stacked vertically, which hold the information for the interpolation of each square. We initialize ZZZ with 0 entries.
    ZZZ = matrix(0,100 * length(squares[,1]),20)
    for(i in 1:length(squares[,1])){
      # we do the triangulation analysis for each square and save it at its index on the ZZZ matrix
      ZZZ[((i - 1) * 100 + 1):(i * 100),] = getTriang(list(x = messSq[1,,i], y = messSq[2,,i]))
    }


    #sq is a vector which assings to each group of 20 in dtSp their respective square index in the squares list
    sq=c()

    for(i in 1:(length(dtSp$t1)/20)){
      sq=c(sq,which(squares[,2]==dtSp$Localidad[i*20] & squares[,3]==dtSp$Exclusion[i*20] & squares[,4]==dtSp$Cuadro[i*20]))
      

    }

    yearmat=as.matrix(dtSp[,10:23])
    11+ncol(yearmat)->firstSp
    ncol(dtSp)->lastSp
    return(list(ZZ=ZZZ,messSq=messSq,sq=sq,yearmat=yearmat))
}

# This function takes a list of two vectors representing the x and y coordinates of the measured subsquares of a square and returns a matrix implicitly encoding a linear interpolation of all 100 z values from the 20 measured ones.
getTriang=function(coord20){
    # 'mkEdg' is a list with the 'x', 'y' and 'z' values of measured subsquares, as well as the necessary information to fill in missing edges of the big square
    mkEdg=makeEdges(coord20$x,coord20$y,rep(0,20))
    coord=list(x=mkEdg$cx,y=mkEdg$cy)
    # Make Delaunay triangulation using 'deldir' package. From it we get the 'sides' of the triangles given by the points they connect.
    delTri=deldir(coord)
    sides=delTri$delsgs
    triangs=triang.list(delTri)
    # 'sideValMat' keeps track, for each triangle side, of whether each of the 100 points is on the 'positive' or 'negative' semiplanes defined by cutting the plane with the line defined by the side. This is better explained in the appendix. 
    sideValMat=matrix(NA,100,length(sides[,1]))
    for(j in 1:length(sides[,1])){
        # We will place each line segment on the x axis with point 1 on (0,0). It's important to keep the order of the points as without this what follows is meaningless...
        vtx=sides[j,]

        # Here we define a translation vector (TT) & rotation matrix (R). TT is actually the first point of each side, and R is the rotatión that makes the side be horizontal
        TT=c(vtx$x1,vtx$y1)
        P2=c(vtx$x2,vtx$y2)
        # Vector from TT to P2 (this is just the side moved to the origin)
        V0=P2-TT
        # The norm of the vector from P1 to P2 (or hypotenuse)
        nrmP1=sqrt(sum(V0**2))
        # R rotational matrix
        R=matrix(c(V0[1]/nrmP1,-V0[2]/nrmP1,+V0[2]/nrmP1,V0[1]/nrmP1),2,2)

        for(i in 1:100){

            # Px takes the values of the i-th subsquares coordinates
            Px=pointToCoord(i)

            # Transformed point
            RTPx=R%*%(Px-TT)

            if(RTPx[2]>0){
                sideValMat[i,j]=1
            }else if(RTPx[2]<0){
                sideValMat[i,j]=-1
            }else{
                sideValMat[i,j]=RTPx[2]
            }
        }
    }

    # Transformed point
    dfPntTrg=data.frame(P=1:100,Triangle=rep(0,100))
    # These three matrices hold the information of which measured subsquares are the vertices of the triangle they belong to. For squares within a triangle whose vertices are not all measured subsquares the idea is to put in the fractions of the measured subsquares that make up the value of the missing subsquare. 
    Z1=matrix(c(0),100,20)
    Z2=matrix(c(0),100,20)
    Z3=matrix(c(0),100,20)
    # 'mtEq' holds the value of the variables used to do the interpolation as explained in the appendix    
    mtEq=matrix(c(0),100,2)

    for(C in 1:length(triangs)){

        t=triangs[[C]]
       
        # Based on the fact that a line won't appear twice on the side list. s1,2,3 are the values we'd expect on the sideVal matrix if the point were towards the inside of a triangle that had that line. s1,s2,s3 define the direction of the side in the triangle with respect to how it appears in sides, and also has the index as its absolute value. 
        #  'max' function is used because 'which' returns 'integer(0)' when there are no matches, and a negative value will be greater for 'max'        
        s1=max(which(sides$ind1==t$ptNum[1] & sides$ind2==t$ptNum[2]),-which(sides$ind1==t$ptNum[2] & sides$ind2==t$ptNum[1]))
        s2=max(which(sides$ind1==t$ptNum[2] & sides$ind2==t$ptNum[3]),-which(sides$ind1==t$ptNum[3] & sides$ind2==t$ptNum[2]))
        s3=max(which(sides$ind1==t$ptNum[3] & sides$ind2==t$ptNum[1]),-which(sides$ind1==t$ptNum[1] & sides$ind2==t$ptNum[3]))

        #Gives a value between 0 and 3 which tells us if the point is inside the triangle C (=3)
        oriSum=abs(rowSums(sideValMat[,abs(c(s1,s2,s3))]*matrix(c(s1/abs(s1),s2/abs(s2),s3/abs(s3)),100,3,byrow=T))) + rowSums(sideValMat[,abs(c(s1,s2,s3))]==matrix(0,100,3))
        #  Ids of points that are inside of 't'        
        inPoints=which(oriSum==3)
        dfPntTrg[inPoints,2]=C

        #Matrices that hold info on which three points determine the triangle inside which each point is in
        Z1[inPoints,]=0
        Z2[inPoints,]=0
        Z3[inPoints,]=0
        
        # If the id of the x point in the triangle is not one of the 20 then the closest n points get value 1/n.
        if(t$ptNum[[1]]<21){
            Z1[inPoints,t$ptNum[[1]]]=1
        }
        else{
            #  If the id of one of the vertices of the triangle is not one of the 20 then the closest n points get value 1/n.
            id="ctl"
            if(coord$x[t$ptNum[[1]]]==1 & coord$y[t$ptNum[[1]]]==10) id="ctr"
            if(coord$x[t$ptNum[[1]]]==10 & coord$y[t$ptNum[[1]]]==1) id="cbl"
            if(coord$x[t$ptNum[[1]]]==10 & coord$y[t$ptNum[[1]]]==10) id="cbr"
            Z1[inPoints,mkEdg[[id]]]=1/length(mkEdg[[id]])
        }
        
        if(t$ptNum[[2]]<21){
            Z2[inPoints,t$ptNum[[2]]]=1
        }
        else{
            id="ctl"
            if(coord$x[t$ptNum[[2]]]==1 & coord$y[t$ptNum[[2]]]==10) id="ctr"
            if(coord$x[t$ptNum[[2]]]==10 & coord$y[t$ptNum[[2]]]==1) id="cbl"
            if(coord$x[t$ptNum[[2]]]==10 & coord$y[t$ptNum[[2]]]==10) id="cbr"
            Z2[inPoints,mkEdg[[id]]]=1/length(mkEdg[[id]])
        }

        if(t$ptNum[[3]]<21){
            Z3[inPoints,t$ptNum[[3]]]=1
        }
        else{
            id="ctl"
            if(coord$x[t$ptNum[[3]]]==1 & coord$y[t$ptNum[[3]]]==10) id="ctr"
            if(coord$x[t$ptNum[[3]]]==10 & coord$y[t$ptNum[[3]]]==1) id="cbl"
            if(coord$x[t$ptNum[[3]]]==10 & coord$y[t$ptNum[[3]]]==10) id="cbr"
            Z3[inPoints,mkEdg[[id]]]=1/length(mkEdg[[id]])
        }

        # Edge points taking into account translation to origin
        TT=as.numeric(t[1,][c(2,3)])
        P2=as.numeric(t[2,][c(2,3)]-TT)
        P3=as.numeric(t[3,][c(2,3)]-TT)
        xy=matrix(pointToCoord(inPoints),length(inPoints),2)
        xy[,1]=xy[,1]-TT[1]
        xy[,2]=xy[,2]-TT[2]
        mtEq[inPoints,]=matrix(c(P2[2]*xy[,1]-P2[1]*xy[,2],-P3[2]*xy[,1]+P3[1]*xy[,2]),length(inPoints),2)/(P2[2]*P3[1]-P2[1]*P3[2])
    }

    finRes=Z1
    for(i in 1:100){
        finRes[i,]=finRes[i,]+(Z2[i,]-Z1[i,])*mtEq[i,2]+(Z3[i,]-Z1[i,])*mtEq[i,1]
    }

    return(finRes)
}

# Adds edges (if missing) of a square to a couple of vectors representing x and y coordinates. Also sets value at missing subsquares as the value of the closest measured subsquares.
makeEdges=function(cx,cy,zz){
    ctl=NaN
    ctr=NaN
    cbl=NaN
    cbr=NaN
    #Top Left
    if(sum(which(cx==1 & cy==1))==0){
        ctl=which(sqrt(cx**2+cy**2)==min(sqrt(cx**2+cy**2)))
        zz=c(zz,mean(zz[ctl]))
        cx=c(cx,1)
        cy=c(cy,1)
    }
    #Top Right
    if(sum(which(cx==1 & cy==10))==0){
        ctr=which((sqrt(cx**2+(10-cy)**2))==min(sqrt(cx**2+(10-cy)**2)))
        zz=c(zz,mean(zz[ctr]))
        cx=c(cx,1)
        cy=c(cy,10)
    }
    #Bottom Left
    if(sum(which(cx==10 & cy==1))==0){
        cbl=which(sqrt((10-cx)**2+cy**2)==min(sqrt((10-cx)**2+cy**2)))
        zz=c(zz,mean(zz[cbl]))
        cx=c(cx,10)
        cy=c(cy,1)
    }
    #Bottom Right
    if(sum(which(cx==10 & cy==10))==0){
        cbr=which(sqrt((10-cx)**2+(10-cy)**2)==min(sqrt((10-cx)**2+(10-cy)**2)))
        zz=c(zz,mean(zz[cbr]))
        cx=c(cx,10)
        cy=c(cy,10)
    }
    return(list(cx=cx,cy=cy,zz=zz,ctl=ctl,ctr=ctr,cbl=cbl,cbr=cbr))
}

pointToCoord = function(i,inv=F){
    if(inv){
        return((i[[2]]-1)*10+(i[[1]]))
    }else{
        return(c(i-10*floor((i-1)/10),floor((i-1)/10)+1))        
    }
}

runCTest=function(dtVeros,paramVeros,model='IXexp3poissinint',ranEf=NULL){
	obj = MakeADFun(dtVeros, paramVeros, DLL=model, random=ranEf)
	obj$hessian = TRUE #????
	opt = nlminb(start = obj$par, obj = obj$fn, gr = obj$gr,control=list(iter.max=1000,eval.max=1000))
	return(opt)

}

runCTestSampled=function(dtVeros,model='IXexp3poissinint',points=10,paramNum=3,ranEf=NULL){
	maxOpt=list() # Here we'll save the best optimization result found
	yrNum=ncol(dtVeros$yearmat)
	spNum=ncol(dtVeros$tx)
	locNum=length(unique(as.numeric(dtSp$Localidad)))
	paramNum=paramNum+yrNum+((!is.null(ranEf))*(locNum+1)) # basic, plus year variable param, plus random effects
	lstValues=c() # A list of the likelyhoods for posterior analysis


	hcPoints=optimumLHS(points,paramNum) # We generate a matrix of points to take as starting points for the optimization process
	hcPoints[,1:yrNum]=(hcPoints[,1:yrNum]-.5)*6 # BB between -3 and 3??
	hcPoints[,yrNum+1]=(hcPoints[,yrNum+1]-.5)*10 # DD between -5 and 5
	hcPoints[,yrNum+2]=(hcPoints[,yrNum+2]*5)-.5 # alfa (values <-0.5  can return nan)
	hcPoints[,yrNum+3]=(hcPoints[,yrNum+3]*4)-2 # cc

	if(!is.null(ranEf)){
		hcPoints[,(paramNum-locNum):(paramNum-1)]=(hcPoints[,(paramNum-locNum):(paramNum-1)]-.5)*10 # Random effects per locality
		hcPoints[,paramNum]=(hcPoints[,paramNum]-.5)*10 # Random variation in RE
	}
	count=1
	crashParams=list()
	for(p in 1:points){
		paramVeros=list(BB=c(hcPoints[p,1:yrNum]),DD=hcPoints[p,yrNum+1],alfa=hcPoints[p,yrNum+2],cc=hcPoints[p,yrNum+3])

		# Random effects
		if(!is.null(ranEf)){
			paramVeros[['ranef']]=hcPoints[p,(paramNum-locNum):(paramNum-1)] # Random effects per locality
			paramVeros[['desvran']]=hcPoints[p,paramNum] # Random variation in RE
		}

		tmpOpt=try(runCTest(dtVeros,paramVeros,model=model,ranEf=ranEf))
		if(class(tmpOpt)=='try-error'){
			crashParams=c(crashParams,paramVeros)
		}else if(length(maxOpt)==0 & tmpOpt$objective>0){
			maxOpt[[count]]=tmpOpt
			lstValues=c(lstValues,tmpOpt$objective)
			count=count+1
		}else #if(isDiferent(tmpOpt,maxOpt) & tmpOpt$objective>0)
		{
			maxOpt[[count]]=tmpOpt
			lstValues=c(lstValues,tmpOpt$objective)
			count=count+1

		}
	}
	return(list(maxOpt=maxOpt,lstValues=lstValues,crashParams=crashParams,points=hcPoints))
}



runBasic=function(dtSp,dtVals,points=100){
	dtVeros=list(t1=dtSp$t1,t2=dtSp$t2,x=dtSp$x,y=dtSp$y,sq=dtVals$sq,yearmat=dtVals$yearmat,ZZ=dtVals$ZZ)
	paramVeros=list(BB=rep(1,dim(dtVals$yearmat)[2]),DD=1,alfa=1,cc=1)
	compile("/home/invitado/Escritorio/Ale/capitulo2/Poisson/IXexp3poissinint.cpp")
	dyn.load(dynlib("IXexp3poissinint"))
	optN=runCTestSampled(dtVeros,points=points)
	return(optN)
}

#read data
dtSp=read.csv("Bouchogibbscor.csv")
#clean data
dtSp=cleanDt(dtSp)
#eliminate zeros from data (estimate of a zero will always be a zero)
dtSp=cleanzero(dtSp) 
dtVals=setupVals(dtSp)
dtSp$Localidad=factor(dtSp$Localidad)

#run model
runBasic(dtSp,dtVals,points=5)->resul
#which run had the maximum likelihood?
resul$maxOpt[[which(resul$lstValues==min(resul$lstValues))]]->maxresul
#maximum likelihood run
maxresul