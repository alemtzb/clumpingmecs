//Alejandra Martínez Blancas & Ian Xul 11/04/21 alemtzb@ciencias.unam.mx
//// Model for second phase parameterization implementation for species with presence/absence data

# include <TMB.hpp>
# include <math.h>
# include <iostream>

using namespace romberg;

template<class Type>
Type objective_function<Type>::operator() ()
{
    //Data
    DATA_MATRIX(tx); //species abundance matrix 
    DATA_MATRIX(txmas); //1 + species abundance matrix
    DATA_VECTOR(interdepth); //vector that specifies which species interactions change with soil depth
    DATA_VECTOR(t1); //species abundance in time t
    DATA_VECTOR(t2); //species abundance in time t +1
    DATA_VECTOR(Dsq); //soil depth of each 10 x 10 cm square
    DATA_VECTOR(Dsq2); //squared soil depth of each 10 x 10 cm square
    DATA_FACTOR(x); //x coordinates of square
    DATA_FACTOR(y); //y coordinates of square
    DATA_FACTOR(sq); //sq holds the correspondence of each set of 20 measured subsquares to the submatrix of ZZ which needs to be used for that square
    DATA_MATRIX(yearmat); //binary matrix for year in the study
    DATA_MATRIX(ZZ); // ZZ is a big matrix that has the information, for each square, to do linear interpolation as a matrix multiplication
    DATA_MATRIX(rainMtrxFull); //dispertion data obtained in phase 1
    DATA_VECTOR(BB); //birth rate parameter estimated in phase 1 which will now use to estimate parameters b and d
    DATA_FACTOR(square); //unique identifier of each quadrat in our data set for the calculation of the random effect

    //Parameters
    PARAMETER(DD); //Death rate
    PARAMETER(sigma); //parameter for calculation of random effects
    PARAMETER_VECTOR(BB2); //parameter b
    PARAMETER_VECTOR(BB3); //parameter c
    PARAMETER_VECTOR(alphas); //parameter g matrix
    PARAMETER_VECTOR(betas); //facilitation parameter matrix
    PARAMETER_VECTOR(alphas2); //parameter h matrix
    PARAMETER_VECTOR(ranef); //random effect

    //declarations
    int numdat=t1.size();
    int lengthpar=alphas.size();
    int numyear=BB2.size();
    vector<Type> immigrants(numdat);
    vector<Type> t2est(numdat);
    vector<Type> surv(numdat);
    vector<Type> emigrants(numdat);
    vector<Type> fac(numdat);
    vector<Type> BBpre(numdat);
    vector<Type> BB1(numyear);
    vector<Type> BBs1(numdat);
    vector<Type> BBs2(numdat);
    vector<Type> BBs3(numdat);
    vector<Type> BBs(numdat);
    vector<Type> numdatcol(numdat);
    vector<Type> profporbb(numdat);


    //bounds
    betas=1/(1+exp(-betas));
    DD=exp(DD)/(1+exp(DD));
    BB2=exp(BB2)/(1+exp(BB2))-0.5;
    BB3=-exp(BB3)/(1+exp(BB3))*5;
    vector<Type> ybetas(lengthpar);
    ybetas << 0, betas;


    //Procedure
    
    //Calculation of birthrate parameters b and c
    BBs2=yearmat*BB2;
    BBs3=yearmat*BB3;
    BBs2=BBs2*Dsq;
    BBs3=BBs3*Dsq2;
    BBpre=exp(BBs2+BBs3);
    for(int i =0;i<numyear;i=i+1) {
        numdatcol=yearmat.col(i);
        profporbb=BBpre*numdatcol;
        Type numdatcolsum=numdatcol.sum();
        Type profporbbsum=profporbb.sum();
        BB1(i)=log(numdatcolsum*BB(i))-log(profporbbsum);
        }
    BBs1=yearmat*BB1;
    BBs=BBs1+BBs2+BBs3;

    //Calculation of competition parameters g and h
    vector <Type> desparams=alphas2*interdepth;
    matrix <Type> desparammat(lengthpar,2);
    desparammat.col(0) = alphas;
    desparammat.col(1) =desparams;
    vector <Type> finalalphas(numdat);
    for(int i=0; i<numdat;i=i+1) {
        vector <Type> depths(2);
        depths << 1,Dsq(i);
        //vector <Type> rowalphas = desparammat*log(depths); //feliz
        vector <Type> rowalphas = desparammat*depths; //feliz
        rowalphas = exp(rowalphas);
        vector <Type> rowtx = tx.row(i);
        vector <Type> rowtxalphas = rowtx * rowalphas;
        finalalphas(i) = rowtxalphas.sum();
    }


    
    //Random effects
    Type REResult = 0; // Initialize the result var
    REResult=-sum(dnorm(ranef,Type(0.0),exp(sigma),true));
    vector<Type> rr(numdat);
    for(int i = 0; i < numdat; i++){
        rr(i)=ranef(square(i));
    }

    //Final calculations
    BBs= BBs + rr;
    BBs=exp(BBs);
    // Individuals that survive
    surv=(1-DD)*t1;
    fac=exp(txmas*ybetas);
    emigrants=BBs*t1/(1.0+finalalphas)*fac;
    

    // This array will hold 100 sheets of 10x10 containing the proportion distributed from each square
    array<Type> sResTot(10,10,100);
    sResTot.fill(0);


    // Now we fill in the sResTot array with the dispersal data obtained in phase 1
    int cont=0;

    matrix<Type> sbRainMtrx(10,10);
    sbRainMtrx.fill(0);

    for( int j = 0; j < 10; j=j+1) {
        for( int i = 0; i < 10; i=i+1){
            sbRainMtrx=rainMtrxFull.block((9-j),(9-i),10,10);
            sResTot.col(cont)=sbRainMtrx.array();
            cont+=1;    
        }
    }


    // The borde matrix will hold data on the necessary correction due to the border effect

    matrix<Type> borde(10,10);
    borde.fill(0);

    // Add up the 100 sheets that make up sResTot
    for(int kk=0; kk<100; kk=kk+1){
        borde=borde + sResTot.col(kk).matrix();
    }

    // Max value of borde is the determining factor, it indicates what the rest should be standartized to 
    Type maxuni=borde.array().maxCoeff();

    // Stand in vectors which will be filled with the values for each square as necessary
    vector<int> cx(20); 
    vector<int> cy(20); 
    vector<Type> zz(20); 
    vector<Type> zzi(100); 
    vector<Type> immiRes(20); 

    // sResTotRot makes it easier to acces the 10 100-long hipercolumns due to limitations in the way arrays are managed in eigen
    array<Type> sResTotRot(100,10,10);
    sResTotRot=sResTot.rotate(1);

    // This for goes through all distinct squares
    for( int i=0; i<(numdat/20); i=i+1){
        cx=x.segment((i*20),20);
        cy=y.segment((i*20),20);
        zz=emigrants.segment((i*20),20);

        //Do the interpolation calculations (refer to synthesis.R) ZZ contains several 20x100 matrices, and is subsetted so only those relevant for the quadrat are considered.
        zzi=ZZ.block((sq(i)-1)*100,0,100,20)*zz.matrix();

        //Calculate the total dispersion within the square given the model (totalDisp) and the total production of successfully dispersed individuals within the square (totalProd)
        Type totalDisp = 0; //number of individuals that remain in the quadrat after dispersal after correcting by borde
        for(int rrow=0;rrow<10;rrow=rrow+1){
            for(int ccol=0;ccol<10;ccol=ccol+1){
                totalDisp=totalDisp+((zzi*sResTotRot.col(rrow).col(ccol)).sum()*(maxuni/borde(rrow,ccol)));
            }
        }
        
        Type totalProd = zzi.sum(); //total number of emigrants. totalProd-totalDisp = number of individual that leave the quadrat.

        // Once we have the 100 z values we sum up the contribution to each square from all others and also correct for borders.
        for(int kk=0; kk<20; kk=kk+1){
            immiRes(kk)=(zzi*sResTotRot.col(cx(kk)-1).col(cy(kk)-1)).sum()*(maxuni/borde(cx(kk)-1,cy(kk)-1));
        }

        immiRes=immiRes+(totalProd-totalDisp)/100; 

        immigrants.segment((i*20),20)=immiRes; 

    }


    // Calculate the estimate (what stays plus what comes from outside)
    t2est=1-(1-surv)*exp(-immigrants);

    Type finalResult = 0; // Initialize the result var
    for(int i = 0; i<numdat; i=i+1){
           finalResult -= log(t2est(i))*t2(i)+log(1-t2est(i))*(1-t2(i));
    }

    return REResult+finalResult;
} 
