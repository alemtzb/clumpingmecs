//Alejandra Martínez Blancas & Ian Xul 11/04/21 alemtzb@ciencias.unam.mx
//// Basic model with modified, more explicit, growth model. Instead of lambdas we have a mortality (D) and a birth (B) rate. 

# include <TMB.hpp>
# include <math.h>
# include <iostream>

using namespace romberg;

// Define the basic structure of the 3-D kernel function which will be integrated
template<class Type>
struct multivariate {
	Type cc;
	Type alfa;               // Parameter in integrand
	multivariate(Type cc_, Type alfa_) // Constructor of integrand
	: cc (cc_),
	alfa (alfa_) {}       // Initializer list
	Type operator()           // Evaluate integrand
	(vector<Type> x){
        Type nrm = 2*M_PI*pow(alfa,2)*exp(lgamma(2/cc))/cc;
		return exp(-(pow(sqrt((x(0)*x(0))+(x(1)*x(1)))/alfa,cc)))/nrm;

	}
};


template<class Type>
Type objective_function<Type>::operator() ()
{
    //Data
    DATA_VECTOR(t1); //species abundance in time t
    DATA_VECTOR(t2); //abundance in time t+1
    DATA_FACTOR(x); //x coordinate of square
    DATA_FACTOR(y); //y coordinate of square
    DATA_FACTOR(sq); // sq holds the correspondence of each set of 20 measured subsquares to the submatrix of ZZ which needs to be used for that square
    DATA_MATRIX(yearmat); //binary matrix for year in the study
    DATA_MATRIX(ZZ); // ZZ is a big matrix that has the information, for each square, to do linear interpolation as a matrix multiplication

    //Parameters
    PARAMETER_VECTOR(BB); //Birth
    PARAMETER(DD); //Death
    PARAMETER(alfa); //dispersal parameter asociated with average dispersion
    PARAMETER(cc); //dispersal parameter associated with kernel shape


    //declarations
    int numdat=t1.size();
    int alphas=0.0;
    vector<Type> immigrants(numdat);
    vector<Type> t2est(numdat);
    vector<Type> surv(numdat);
    vector<Type> emigrants(numdat);
    vector<Type> BBs(numdat);

    //bounds
    alfa=exp(alfa);
    cc=(exp(cc)/(1+exp(cc)))*2;
    DD=exp(DD)/(1+exp(DD));
    BBs=exp(yearmat*BB);



    //procedure
    // Individuals that survive
    surv=(1-DD)*t1;
    // Individuals that succesfully disperse
    emigrants=BBs*t1/(1.0+t1*alphas);


    // This array will hold 100 sheets of 10x10 containing the proportion distributed from each square
	array<Type> sResTot(10,10,100);
	sResTot.fill(0);

    // This matrix will be filled with the top half of the dispersion kernel from square 1,1. The other half is just the reflection
	matrix<Type> rainMtrxHalf(10,10);
	rainMtrxHalf.fill(0);

    // Instatiate the multivariate function 'f' and then use it to do integration later on (we only need to instatiate it once)
    multivariate<Type> f(cc,alfa);

	for( int i = 0; i < 10; i=i+1) {
    	for( int j = i; j < 10; j=j+1){
            // We create a couple of vectors which define the square of integration for the dispersion kernel
    		vector<Type> ss(2);
    		vector<Type> tt(2);
    		ss(0)=j;
    		ss(1)=i;
    		tt(0)=j+1;
    		tt(1)=i+1;

            // Integrate
			Type shad = romberg::integrate(f, ss, tt);

    		rainMtrxHalf(i,j)=shad;
    	}
    }
    
	// We set the lower triangular part of rainMtrx equal to the upper. Due to data types it is a little less straight-forward than in R

    matrix<Type> rainMtrx(10,10);
    matrix<Type> rainMtrxHalf2(10,10);

    rainMtrxHalf2=rainMtrxHalf.transpose();
    rainMtrxHalf2.diagonal().fill(0);
    rainMtrx=rainMtrxHalf+rainMtrxHalf2;

    // Fill in the rest of the 19x19 kernel defining matrix by reflecting and pasting rainMtrx 

    matrix<Type> rainMtrxFull(19,19);
    rainMtrxFull.fill(0);

    rainMtrxFull.bottomRightCorner(10,10)=rainMtrx;
    rainMtrxFull.topRightCorner(9,10)=rainMtrx.colwise().reverse().topRows(9);
    rainMtrxFull.topLeftCorner(10,9)=rainMtrx.reverse().leftCols(9);
    rainMtrxFull.bottomLeftCorner(9,9)=rainMtrx.rowwise().reverse().bottomLeftCorner(9,9);


    // Now we fill in the sResTot array with each 'sheet' selected as a 10x10 matrix from the 19x19 mt centered on the interest square

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

    // sResTotRot makes it easier to acces the 10 100-long hipercolumns due to limitations in the way arrays are managed in eigen
    array<Type> sResTotRot(100,10,10);
    sResTotRot=sResTot.rotate(1);

    // Add up the 100 sheets that make up sResTot to obtain the borde matrix
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

    // This for goes through all distinct squares
    for( int i=0; i<(numdat/20); i=i+1){
    	cx=x.segment((i*20),20);
    	cy=y.segment((i*20),20);
    	zz=emigrants.segment((i*20),20);

        // Do the interpolation calculations (refer to the appendix) ZZ contains 20x100 matrices which when multiplied by zz give the interpolated z values
    	zzi=ZZ.block((sq(i)-1)*100,0,100,20)*zz.matrix();

        //(IX) Calculate the total dispersion within the square given the model (totalDisp) and the total production of successfully dispersed individuals within the square (totalProd)    
        Type totalDisp = 0;
        for(int rr=0;rr<10;rr=rr+1){
            for(int cc=0;cc<10;cc=cc+1){
                totalDisp=totalDisp+((zzi*sResTotRot.col(rr).col(cc)).sum()*(maxuni/borde(rr,cc)));
            }
        }
        totalDisp=totalDisp;

        Type totalProd = zzi.sum();


        // Once we have the 100 z values we sum up the contribution to each square from all others and also correct for borders.
        for(int kk=0; kk<20; kk=kk+1){
            immiRes(kk)=(zzi*sResTotRot.col(cx(kk)-1).col(cy(kk)-1)).sum()*(maxuni/borde(cx(kk)-1,cy(kk)-1));
            immiRes(kk)=immiRes(kk)+(totalProd-totalDisp)/100;
        }


    	immigrants.segment((i*20),20)=immiRes;

    }

    // Calculate the probability that the square remains ocupied (what stays plus what comes from outside) or is conlonized assuming a poisson distribution
    t2est=1-(1-surv)*exp(-immigrants);

    Type finalResult = 0; // Initialize the result var


    // Best I could think of but would be better just to filter out. Sums up prob. only for datum with estimates different from 0 (otherwise it doesn't make sense)

    for(int i = 0; i<numdat; i=i+1){
    	   finalResult -=log(t2est(i))*t2(i)+log(1-t2est(i))*(1-t2(i));
    }
    

    return finalResult;
    //Format for reporting variables
}
