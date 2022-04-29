# include <TMB.hpp>
# include <math.h>
# include <iostream>

using namespace romberg;

template<class Type>
Type objective_function<Type>::operator() ()
{
    //Data
    DATA_MATRIX(tx);
    DATA_MATRIX(txmas);
    DATA_VECTOR(interdepth);
    DATA_VECTOR(t1);
    DATA_VECTOR(t2);
    DATA_VECTOR(Dsq);
    DATA_VECTOR(Dsq2);
    DATA_FACTOR(x);
    DATA_FACTOR(y);
    DATA_FACTOR(sq);
    DATA_MATRIX(yearmat);
    DATA_MATRIX(ZZ);
    DATA_MATRIX(rainMtrxFull);
    DATA_VECTOR(BB);
    DATA_FACTOR(square);

    //Parameters
    PARAMETER(DD);
    PARAMETER(theta);
    PARAMETER(sigma); 
    PARAMETER_VECTOR(BB2);
    PARAMETER_VECTOR(BB3);
    PARAMETER_VECTOR(alphas);
    PARAMETER_VECTOR(betas);
    PARAMETER_VECTOR(alphas2);    
    //PARAMETER_VECTOR(betas2);
    PARAMETER_VECTOR(ranef);

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
    //alphas=exp(alphas);
    betas=1/(1+exp(-betas));
    theta=exp(theta);
    DD=exp(DD)/(1+exp(DD));
    BB2=exp(BB2)/(1+exp(BB2))-0.5;
    BB3=-exp(BB3)/(1+exp(BB3))*5;
    vector<Type> ybetas(lengthpar);
    ybetas << 0, betas;
    
    //vector<Type> ybetas2(lengthpar);
    //ybetas2 << 0, betas2;

    //Procedure
    
    //Birth
    
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

    //alphas
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

    //Final calculation
    BBs= BBs + rr;
    BBs=exp(BBs);

    surv=(1-DD)*t1;
    fac=exp(txmas*ybetas);
    emigrants=BBs*t1/(1.0+finalalphas)*fac;
    //emigrants=BBs*t1/(1.0+finalalphas);
    
    //vector<Type> fac_head = fac.head(6);
    //Rcout<<fac_head<<std::endl;

// This array will hold 100 sheets of 10x10 containing the proportion distributed from each square
    array<Type> sResTot(10,10,100);
    sResTot.fill(0);

    // Now we fill each sheet in the sResTot by cutting out a 10x10 section of rainMtrxFull so that the source square is placed in the respective position.

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

    // Max value of borde is the reference value, it indicates what the rest should be standardized to 
    Type maxuni=borde.array().maxCoeff();

    //3. Estimate dispersal for each quadrat
    // Declare vectors which will be filled with the values for each square as necessary
    vector<int> cx(20); // x coordinates of the 20 observed squares in the quadrat
    vector<int> cy(20); // y coordinates of the 20 observed squares in the quadrat
    vector<Type> zz(20); // number of emigrants from the 20 observed squares in the quadrat
    vector<Type> zzi(100); //number of emigrants from all the 100 squares in the quadrat
    vector<Type> immiRes(20); //number of immigrants into the 20 observed squares in the quadrat

    // sResTotRot makes it easier to acces the 10 100-long hipercolumns due to limitations in the way arrays are managed in eigen
    array<Type> sResTotRot(100,10,10);
    sResTotRot=sResTot.rotate(1);

    // This for goes through all squares
    for( int i=0; i<(numdat/20); i=i+1){
        //Select the data for each square
        cx=x.segment((i*20),20);
        cy=y.segment((i*20),20);
        zz=emigrants.segment((i*20),20);

        //Do the interpolation calculations (refer to synthesis.R) ZZ contains several 20x100 matrices, and is subsetted so only those relevant for the quadrat are considered.
        //These are multiplied by zz give the interpolated number of emigrants from every square in the matrix.
        zzi=ZZ.block((sq(i)-1)*100,0,100,20)*zz.matrix();

        Type totalDisp = 0; //number of individuals that remain in the quadrat after dispersal after correcting by borde
        for(int rrow=0;rrow<10;rrow=rrow+1){
            for(int ccol=0;ccol<10;ccol=ccol+1){
                totalDisp=totalDisp+((zzi*sResTotRot.col(rrow).col(ccol)).sum()*(maxuni/borde(rrow,ccol)));
            }
        }
        
        Type totalProd = zzi.sum(); //total number of emigrants. totalProd-totalDisp = number of individual that leave the quadrat.

        // Do the same as before (sum up the contribution to each square from all others and also correct for borders), but only for the 20 observed squares values.
        for(int kk=0; kk<20; kk=kk+1){
            immiRes(kk)=(zzi*sResTotRot.col(cx(kk)-1).col(cy(kk)-1)).sum()*(maxuni/borde(cx(kk)-1,cy(kk)-1));
        }

        // Assuming that the number of seeds that immigrate to the quadrat is same number that emigrate out of it, and redistributing them equally among all the squares:
        immiRes=immiRes+(totalProd-totalDisp)/100; //immiRes has the number of seeds that arrive to each of the 20 squares in the quadrat

        immigrants.segment((i*20),20)=immiRes; //fill up immigrants (which contains data from all quadrats)

    }


    // Calculate the estimate (what stays plus what comes from outside)
    //t2est=surv+immigrants;
    t2est=surv+immigrants;

    // Best I could think of but would be better just to filter out. Sums up prob. only for datum with estimates different from 0 (otherwise it doesn't make sense)
    Type finalResult = 0; // Initialize the result var
    for(int i = 0; i<numdat; i=i+1){
        //if(t2est(i)>0){
           finalResult -= dnbinom(t2(i), theta, theta/(theta+t2est(i)), 1);
        //}
    }
    //Rcout<<(theta/(theta+t2est(1)))<<t2est(1)<<dnbinom(t2(1), theta, theta/(theta+t2est(1)), 1)<<std::endl;
    //Rcout<<REResult+finalResult<<std::endl;
    return REResult+finalResult;


    //Format for reporting variables
    //Rcout<<t2(i)<<' '<<t2est(i)<<' '<<dnbinom(t2(i), theta, theta/(theta+t2est(i)), 1)<<std::endl;
}
