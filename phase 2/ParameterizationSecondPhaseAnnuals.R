#Alejandra Martínez Blancas & Ian Xul 11/04/21 alemtzb@ciencias.unam.mx
require(TMB)
require(lhs)
require(deldir)
require(cubature)

#load parameters calculated en phase 1
read.csv("clumpingmecs/phase 2/disp.csv")->disp #seed movement parameters
read.csv("clumpingmecs/phase 2/Phase1BR14y.csv",header=F)->BBinp #Birth rate parameters obtained in Phase 1
BBinp=data.frame(BBinp,row.names=T)

cleanDt=function(dtSp,lagNum=0){
    #Remove these combinations of site/year:
    ## If I wrote the names well there should be no problems
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

#Eliminates quadrats that are completely void from plants
cleanzero=function(dtSp){
    numcuad=nrow(dtSp)/20
    elim=c()
    for(i in 1:numcuad){
            if(sum(dtSp$t1[((i-1)*20+1):(i*20)])==0) elim=c(elim,((i-1)*20+1):(i*20))
    }
    return(dtSp[-elim,])
}

setupVals = function(dtSp,disp=disp,BBinp=BBinp,especie,lagNum=0){# dtSp=read.csv(file.choose())
    #squares list will hold the different combinations of factors that determine the squares
    squares=list()
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

    # ZZZ holds n 100 long matrices which hold the information for the triangulation of each 
    ZZZ=matrix(0,100*length(squares[,1]),20)
    for(i in 1:length(squares[,1])){
      # we do the triangulation for each case and hold it temporarily
      ZZZ[((i-1)*100+1):(i*100),]=getTriang(list(x=messSq[1,,i],y=messSq[2,,i]))
    }

    #sq is a vector which assings to each group of 20 in dtSp their respective square index in the squares list
    sq=c()
    #sqMap is a list with the number, for each square, of its primitives (same square one year before) location
    sqMap=list()
    for(lg in (0:lagNum)[-1]){
        # All values are greater than one if the primitive exists. If primitive number is zero it is because the values there were of zero, so it was removed.
        sqMap[[paste('tm',as.character(lg-1),sep='')]]=rep(0,length(dtSp$t1)/20)
    }
    for(i in 1:(length(dtSp$t1)/20)){
      sq=c(sq,which(squares[,2]==dtSp$Localidad[i*20] & squares[,3]==dtSp$Exclusion[i*20] & squares[,4]==dtSp$Cuadro[i*20]))
      for(lg in (0:lagNum)[-1]){
        for(wch in which(dtSp$Localidad==dtSp$Localidad[1+(i-1)*20] & dtSp$Exclusion==dtSp$Exclusion[1+(i-1)*20] & dtSp$Cuadro==dtSp$Cuadro[1+(i-1)*20] & dtSp$year==dtSp$year[1+(i-1)*20]-1)){
            sqMap[[paste('tm',as.character(lg-1),sep='')]][i]=1+(wch-1)/20
            break
        }
      }
    }
    
    yearmat=as.matrix(dtSp[,10:23])
    11+ncol(yearmat)->firstSp
    ncol(dtSp)->lastSp
    tx=as.matrix(dtSp[,c(8,firstSp:lastSp)])
    (ncol(tx)-1)/2->spnum
    tx[,spnum+2:(spnum+1)]->interdepth
    as.vector(interdepth[1,])->interdepth
    interdepth=c(0,interdepth)
    tx[,1:(spnum+1)]->tx
    
    which(disp$X==especie)->sp
    disp[sp,2]->alfa
    disp[sp,3]->cc
    rainMtrx=generateRainMtrx(alfa,cc)
    
    which(rownames(BBinp)==especie)->sp
    BBinp[sp,]->BB
    as.matrix(BB)[1,]->BB
    
    matrix(nrow=nrow(tx),ncol=ncol(tx))->DsqMat
    dtSp$Profundidad->DsqMat[,1:ncol(DsqMat)]
    
    square=as.numeric(as.factor(paste(dtSp$Localidad,dtSp$Exclusion,dtSp$Cuadro))) 
    return(list(ZZ=ZZZ,messSq=messSq,sq=sq,yearmat=yearmat,sqMap=sqMap,tx=tx,interdepth=interdepth,rain=rainMtrx,square=square-1,BB=BB))
}


# This function takes a list of two vectors representing the x and y coordinates of the measured subsquares of a square and returns a matrix implicitly encoding a linear interpolation of all 100 z values from the 20 measured ones.
getTriang=function(coord20){
    mkEdg=makeEdges(coord20$x,coord20$y,rep(0,20))
    coord=list(x=mkEdg$cx,y=mkEdg$cy)
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

    # 'dfPntTrg' will hold, for each of the 100 subsquares, the ID of one of the triangles it belongs to (if it is on the side of a triangle it will belong to two, but it doesn't matter which one you choose).
    dfPntTrg=data.frame(P=1:100,Triangle=rep(0,100))
    # These three matrices hold the information of which measured subsquares are the vertices of the triangle they belong to. For squares within a triangle whose vertices are not all measured subsquares the idea is to put in the fractions of the measured subsquares that make up the value of the missing subsquare. 
    Z1=matrix(c(0),100,20)
    Z2=matrix(c(0),100,20)
    Z3=matrix(c(0),100,20)
    mtEq=matrix(c(0),100,2)

    for(C in 1:length(triangs)){

        t=triangs[[C]]
       
        # Based on the fact that a line won't appear twice on the side list. s1,2,3 are the values we'd expect on the sideVal matrix if the point were towards the inside of a triangle that had that line. s1,s2,s3 define the direction of the side in the triangle with respect to how it appears in sides, and also has the index as its absolute value. 
        s1=max(which(sides$ind1==t$ptNum[1] & sides$ind2==t$ptNum[2]),-which(sides$ind1==t$ptNum[2] & sides$ind2==t$ptNum[1]))
        s2=max(which(sides$ind1==t$ptNum[2] & sides$ind2==t$ptNum[3]),-which(sides$ind1==t$ptNum[3] & sides$ind2==t$ptNum[2]))
        s3=max(which(sides$ind1==t$ptNum[3] & sides$ind2==t$ptNum[1]),-which(sides$ind1==t$ptNum[1] & sides$ind2==t$ptNum[3]))

        #Gives a value between 0 and 3 which tells us if the point is inside the triangle C (=3)
        oriSum=abs(rowSums(sideValMat[,abs(c(s1,s2,s3))]*matrix(c(s1/abs(s1),s2/abs(s2),s3/abs(s3)),100,3,byrow=T))) + rowSums(sideValMat[,abs(c(s1,s2,s3))]==matrix(0,100,3))
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
        #return((i[[1]]-1)*10+(i[[2]]))
        return((i[[2]]-1)*10+(i[[1]]))
    }else{
        #return(c(floor((i-1)/10)+1,i-10*floor((i-1)/10)))        
        return(c(i-10*floor((i-1)/10),floor((i-1)/10)+1))        
    }
}

#uses dispersion parameters obtained in phase 1 to compute seed movement
generateRainMtrx = function(alfa,cc){
    rainMtrx=matrix(0,10,10)
    tar<-c(0,0)
    #Double for sifts only through the upper triangle of sRes. Produces upper triangular probability matrix.
    for(i in 1:10){
            for(j in i:10){
                    ss <- c(j-1,i-1)

                    pp <- abs(ss-tar)
                    # hcubature integrates the kernel for a given hipercube(square in this case). I compressed the whole kernel function into one line.
                    shad=hcubature(function(vv,alfa,cc){
                        rr <- sqrt((vv[1]**2)+(vv[2]**2));
                        nrm <- 2*pi*(alfa**2)*gamma(2/cc)/cc;
                        return(exp(-((rr/alfa)**cc))/nrm);
                    },(pp),(pp+1),alfa=alfa,cc=cc,maxEval = 100 )

                    # print(paste("Shadow from ",pp[1]," ",pp[2]," to ",pp[1]+1," ",pp[2]+1," : ",shad$integral))

                    rainMtrx[i,j] <- abs(shad$integral)
            }
    }
    #rainMtrx2 is used to fill in lower triangle of matrix.
    rainMtrx2=rainMtrx
    diag(rainMtrx2)=0
    rainMtrx=rainMtrx+t(rainMtrx2)

    rainMtrx=cbind(rainMtrx[,10:2],rainMtrx)
    rainMtrx=rbind(rainMtrx[10:2,],rainMtrx)
}

runCTest=function(dtVeros,paramVeros,model='ParameterizationSecondPhaseAnnuals'){
	obj = MakeADFun(dtVeros, paramVeros, DLL=model,random="ranef") #ranEf=ranEf con efectos aleatorios
	obj$hessian = TRUE #????
	opt = nlminb(start = obj$par, obj = obj$fn, gr = obj$gr,control=list(iter.max=10000,eval.max=10000))
	return(opt)

}

runCTestSampled=function(dtVeros,model='ParameterizationSecondPhaseAnnuals',points=10,paramNum=2,interdepth=dtVals$interdepth){
    maxOpt=list() # Here we'll save the best optimization result found
    yrNum=ncol(dtVeros$yearmat)
    spNum=ncol(dtVeros$tx)
    sqNum=length(unique(dtVals$square))
    alphas2=which(interdepth==1)
    nalphas2=length(alphas2)
    if(nalphas2==0){ #if the competition with none of the species changes with soil depth
        paramNum=paramNum+((2*yrNum)+(2*spNum)-1) #+sqNum # basic, plus year variable param, plus random effects
    lstValues=c() # A list of the likelyhoods for posterior analysis

    hcPoints=geneticLHS(points,paramNum) # We generate a matrix of points to take as starting points for the optimization process. The points are bounded.
    hcPoints[,1]=(hcPoints[,1]*2)-1 # Theta 
    hcPoints[,2]=(hcPoints[,2]*2)-1 #Sigma
    hcPoints[,3:(2+yrNum)]=(hcPoints[,3:(2+yrNum)]*2)-1 parameter b
    hcPoints[,(3+yrNum):(2+yrNum*2)]=(hcPoints[,(3+yrNum):(2+yrNum*2)]*2)-1 #parameter c
    hcPoints[,(3+yrNum*2):(2+yrNum*2+spNum)]=(hcPoints[,(3+yrNum*2):(2+yrNum*2+spNum)]*2)-1 # parameter g
    hcPoints[,(3+yrNum*2+spNum):(1+yrNum*2+spNum*2)]=(hcPoints[,(3+yrNum*2+spNum):(1+yrNum*2+spNum*2)]*2)-1 #facilitation parameter
    matre=matrix(runif(points*sqNum)*2-1,nrow=points,ncol=length(unique(dtVals$square)))

    count=1
    crashParams=list()
    for(p in 1:points){
        paramVeros=list(theta=hcPoints[p,1],sigma=hcPoints[p,3],BB2=c(hcPoints[p,3:(2+yrNum)]),BB3=c(hcPoints[p,(3+yrNum):(2+yrNum*2)]),alphas=c(hcPoints[p,(3+yrNum*2):(2+yrNum*2+spNum)]),betas=c(hcPoints[p,(3+yrNum*2+spNum):(1+yrNum*2+spNum*2)]),alphas2=rep(0,ncol(dtVals$tx)),ranef=matre[p,])


        tmpOpt=try(runCTest(dtVeros,paramVeros,model=model))
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
    return(list(maxOpt=maxOpt,lstValues=lstValues,crashParams=crashParams,points=hcPoints,RE=matre))

    }
    else { #if competition with at least one species changes with soil depth
        paramNum=paramNum+((2*yrNum)+(2*spNum)+nalphas2-1) #+sqNum # basic, plus year variable param, plus random effects
    lstValues=c() # A list of the likelyhoods for posterior analysis

    hcPoints=geneticLHS(points,paramNum) # We generate a matrix of points to take as starting points for the optimization process
    hcPoints[,1]=(hcPoints[,1]*2)-1 # Theta 
    hcPoints[,2]=(hcPoints[,2]*2)-1 #Sigma
    hcPoints[,3:(2+yrNum)]=(hcPoints[,3:(2+yrNum)]*2)-1 #parameter b
    hcPoints[,(3+yrNum):(2+yrNum*2)]=(hcPoints[,(3+yrNum):(2+yrNum*2)]*2)-1 #parameter c
    hcPoints[,(3+yrNum*2):(2+yrNum*2+spNum)]=(hcPoints[,(3+yrNum*2):(2+yrNum*2+spNum)]*2)-1 # parameter g 
    hcPoints[,(3+yrNum*2+spNum):(1+yrNum*2+spNum*2)]=(hcPoints[,(3+yrNum*2+spNum):(1+yrNum*2+spNum*2)]*2)-1 #facilitation parameters
    hcPoints[,(2+yrNum*2+spNum*2):(1+(yrNum*2)+(spNum*2)+nalphas2)]=(hcPoints[,(2+yrNum*2+spNum*2):(1+(yrNum*2)+(spNum*2)+nalphas2)]*2)-1 #parameter h
    matre=matrix(runif(points*sqNum)*2-1,nrow=points,ncol=length(unique(dtVals$square)))

    count=1
    crashParams=list()
    for(p in 1:points){
        hcPoints[p,(2+yrNum*2+spNum*2):(1+(yrNum*2)+(spNum*2)+nalphas2)]->prealphas2
        for(n in 1:nalphas2){
            interdepth[alphas2[n]]=prealphas2[n]
        }
        paramVeros=list(theta=hcPoints[p,1],sigma=hcPoints[p,2],BB2=c(hcPoints[p,3:(2+yrNum)]),BB3=c(hcPoints[p,(3+yrNum):(2+yrNum*2)]),alphas=c(hcPoints[p,(3+yrNum*2):(2+yrNum*2+spNum)]),betas=c(hcPoints[p,(3+yrNum*2+spNum):(1+yrNum*2+spNum*2)]),alphas2=c(interdepth),ranef=matre[p,])


        tmpOpt=try(runCTest(dtVeros,paramVeros,model=model))
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
    return(list(maxOpt=maxOpt,lstValues=lstValues,crashParams=crashParams,points=hcPoints,RE=matre))
}
}

runBasic=function(dtSp,dtVals,points=100){
	dtVeros=list(tx=dtVals$tx,txmas=log(dtVals$tx+1),interdepth=dtVals$interdepth,t1=dtSp$t1,t2=dtSp$t2,Dsq=dtSp$Profundidad,Dsq2=(dtSp$Profundidad^2)/1000,x=dtSp$x,y=dtSp$y,sq=dtVals$sq,yearmat=dtVals$yearmat,ZZ=dtVals$ZZ,rainMtrxFull=dtVals$rain,BB=dtVals$BB,square=dtVals$square,DD=1)
paramVeros=list(theta=0.9358384,sigma=1,BB2=rep(1,dim(dtVals$yearmat)[2]),BB3=rep(1,dim(dtVals$yearmat)[2]),alphas=c(0.2254611,rep(1,(ncol(dtVals$tx)-1))),betas=rep(1,ncol(dtVals$tx)-1),alphas2=rep(1,ncol(dtVals$tx)),ranef=rep(0,length(unique(dtVals$square))))#sigma=0,ranef=rep(0,length(unique(dtVals$square))))
    compile("/clumpingmecs/phase 2/ParameterizationSecondPhaseAnnuals.cpp")
    dyn.load(dynlib("ParameterizationSecondPhaseAnnuals"))
	optN=runCTestSampled(dtVeros,points=points)
	return(optN)
}

#read data
dtSp=read.csv("HetpinintMegasinMuhphaBultensinBouchoprof.csv") #specify the csv file for one of the focal species
#clean data
dtSp=cleanDt(dtSp)
#eliminate zeros from data (estimate of a zero will always be a zero)
dtSp=cleanzero(dtSp) 
dtVals=setupVals(dtSp,disp=disp,BBinp=BBinp,"Hetpin")#the third parameter of this function must be the first three letters of the generic name and the first thre letters of the specific epithet of the species. For example Micrhochlo kunthii would be "Mickun". This is to associate the dispersal parameters calculated in the first phase with the parameters being esttimated in this phase for each species 

#run model
runBasic(dtSp,dtVals,points=500)->resul
#save results
saveRDS(resul,"500hetpinsinMuhphaBultensinBouchoprof")
#which run had the maximum liklihood
resul$maxOpt[[which(resul$lstValues==min(resul$lstValues))]]->maxresul
#maximum liklihhood run
maxresul

