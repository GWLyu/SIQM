//
//  BHNSNumericalEvolver.c
//  
//
//  Created by Michael Lahaye on 2022-04-02.
//

#include "BHNSNumericalEvolver.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "BHNSVectors.h"
#include "BHNSWaveform.h"
#include <complex.h>



//Function to compute Omega1
//Inputs:
//      -the quadrupole parameter of the body
//      -the mass of the body
//      -the current value of the PN parameter y
//      -the mass normalized spin vector of the body
//      -the mass normalized spin vector of the other body
//      -the normalized orbital angular momentum vector Lhat
//      -location to return the result
//Outputs:
//      -none
void CalcOmega(long double kappa , long double m , long double y , long double s1[3] , long double s2[3] , long double Lhat[3] , long double returnVector[3] )
{
    long double temp1[3];
    long double temp2[3];
    
    CrossProduct( Lhat , s1 , temp1 );
    long double coeff = ( 0.5*m + 1.5*( 1 - y*( kappa*DotProduct( s1 , Lhat ) + DotProduct( s2 , Lhat ) ) ) );
    ScalarMultiply( coeff , temp1 , temp1 );
    
    CrossProduct( s2 , s1 , temp2 );
    ScalarMultiply( 0.5*y , temp2 , temp2 );
    
    SumVectors( temp1 , temp2 , returnVector );
}



//Function to compute the vector derivatives of s1, s2 and Lhat simultaneously
//Inputs:
//      -the quadrupole parameter of body 1
//      -the quadrupole parameter of body 2
//      -the mass of body 1
//      -the mass of body 2
//      -the current value of the PN parameter y
//      -the mass normalized spin vector s1
//      -the mass normalized spin vector s2
//      -the normalized orbital angular momentum vector Lhat
//      -location to return the result for the derivative of s1
//      -location to return the result for the derivative of s2
//      -location to return the result for the derivative of Lhat
//Outputs:
//      -none
void CalcAngularMomentaDerivs( long double M , long double kappa1 , long double kappa2 , long double m1 , long double m2 , long double y , long double s1[3] , long double s2[3] , long double Lhat[3] , long double returnVectors1[3] , long double returnVectors2[3] , long double returnVectorLhat[3] )
{
    long double temp1[3];
    long double temp2[3];
    
    CalcOmega( kappa1 , m1 , y , s1 , s2 , Lhat , temp1 );
    CalcOmega( kappa2 , m2 , y , s2 , s1 , Lhat , temp2 );
    
    ScalarMultiply( m2*powl( y , 5.0 )/M , temp1 , returnVectors1 );
    ScalarMultiply( m1*powl( y , 5.0 )/M , temp2 , returnVectors2 );
    SumVectors( temp1 , temp2 , returnVectorLhat );
    ScalarMultiply( (-1.0)*powl( y , 6.0 )/M , returnVectorLhat , returnVectorLhat );
}



//Function to compute the derivative of y
long double Calcdydt( long double y , long double eta )
{
    return (32.0*eta/5.0)*powl( y , 9.0);
}


//function to calculate the derivative of phiz
//Inputs:
//      -the mass of body 1
//      -the mass of body 2
//      -the mass normalized spin vector s1
//      -the mass normalized spin vector s2
//      -the normalized orbital angular momentum vector Lhat
//Outputs:
//      -the derivative of phiz based on these
void CalcdPhizdtAnddZetadt(  long double M , long double m1 , long double m2 , long double eta , long double y , long double s1[3] , long double s2[3] , long double Lhat[3] , long double dLhat[3] , long double returnVec[2] )
{
    long double J[3] , cos , temp;
    J[0] = eta*Lhat[0]/y + m1*s1[0] + m2*s2[0];
    J[1] = eta*Lhat[1]/y + m1*s1[1] + m2*s2[1];
    J[2] = eta*Lhat[2]/y + m1*s1[2] + m2*s2[2];
    
    Normalize(J);
    cos = DotProduct( J , Lhat );
    
    CrossProduct( J , Lhat , J );
    
    returnVec[0] = M*DotProduct( dLhat , J )/( 1 - cos*cos )/M;
    returnVec[1] = -cos*returnVec[0];
    
}



//function to calculate deltachi for plotting purposes
//Inputs:
//      -the mass normalized spin vector s1
//      -the mass normalized spin vector s2
//      -the normalized orbital angular momentum vector Lhat
//Outputs:
//      -deltachi
long double CalcDeltachi( long double s1[3] , long double s2[3] , long double Lhat[3] )
{

    return DotProduct( s1 , Lhat ) - DotProduct( s2 , Lhat );
    
}



//function to calculate chieff for plotting purposes
//Inputs:
//      -the mass normalized spin vector s1
//      -the mass normalized spin vector s2
//      -the normalized orbital angular momentum vector Lhat
//Outputs:
//      -chieff
long double CalcChieff( long double s1[3] , long double s2[3] , long double Lhat[3] )
{

    return DotProduct( s1 , Lhat ) + DotProduct( s2 , Lhat );
    
}



//function to calculate J for plotting purposes
//Inputs:
//      -the mass normalized spin vector s1
//      -the mass normalized spin vector s2
//      -the normalized orbital angular momentum vector Lhat
//      -the mass of body 1
//      -the mass of body 2
//      -the current value of the PN parameter y
//      -the symmetric mass ratio
//Outputs:
//      -the magnitude of J
long double CalcJ( long double s1[3] , long double s2[3] , long double Lhat[3] , long double m1 , long double m2 , long double y , long double eta )
{

    long double J[3];
    J[0] = eta*Lhat[0]/y + m1*s1[0] + m2*s2[0];
    J[1] = eta*Lhat[1]/y + m1*s1[1] + m2*s2[1];
    J[2] = eta*Lhat[2]/y + m1*s1[2] + m2*s2[2];
    return powl( J[0]*J[0] + J[1]*J[1] + J[2]*J[2] , 0.5 );
    
}



//Function to compute the vectors to add to s1, s2, Lhat, and the scalars to add to y and phiz according to an RK4 evolution of the derivatives.
//Inputs:
//      -the quadrupole parameter of body 1
//      -the quadrupole parameter of body 2
//      -the mass of body 1
//      -the mass of body 2
//      -the current value of the PN parameter y
//      -the mass normalized spin vector s1
//      -the mass normalized spin vector s2
//      -the normalized orbital angular momentum vector Lhat
//      -the time difference to the next point
//      -array of arrays to return the things to be added to s1, s2, Lhat respectively
//      -array to return the things to be added to y and phiz respectively
void RK4evolver( long double M , long double kappa1 , long double kappa2 , long double m1 , long double m2 , long double y , long double s1[3] , long double s2[3] , long double Lhat[3] , long double deltat , long double angularDerivs[3][3] , long double scalarDerivs[3] )
{
    long double k1vec[3][3] , k2vec[3][3] , k3vec[3][3], k4vec[4][4] , k1scal[3] , k2scal[3] , k3scal[3] , k4scal[3];
    long double temps1[3] , temps2[3], tempLhat[3] , tempy , tempphizandzeta[2];
    
    //calculate k1
    CalcAngularMomentaDerivs( M , kappa1 , kappa2 , m1 , m2 , y , s1 , s2 , Lhat , k1vec[0] , k1vec[1] , k1vec[2] );
    k1scal[0] = Calcdydt( y , m1*m2 )/M;
    CalcdPhizdtAnddZetadt( M , m1 , m2 , m1*m2 , y , s1 , s2 , Lhat , k1vec[2] , tempphizandzeta );
    k1scal[1] = tempphizandzeta[0];
    k1scal[2] = tempphizandzeta[1];
    
    //calculate k2
    ScalarMultiply( deltat*0.5 , k1vec[0] , temps1 );
    SumVectors( temps1 , s1 , temps1 );
    ScalarMultiply( deltat*0.5 , k1vec[1] , temps2 );
    SumVectors( temps2 , s2 , temps2 );
    ScalarMultiply( deltat*0.5 , k1vec[2] , tempLhat );
    SumVectors( tempLhat , Lhat , tempLhat );
    tempy = y + 0.5*deltat*k1scal[0];
    CalcAngularMomentaDerivs( M , kappa1 , kappa2 , m1 , m2 , tempy , temps1 , temps2 , tempLhat , k2vec[0] , k2vec[1] , k2vec[2] );
    k2scal[0] = Calcdydt( tempy , m1*m2 )/M;
    CalcdPhizdtAnddZetadt( M , m1 , m2 , m1*m2 , tempy , temps1 , temps2 , tempLhat , k2vec[2] , tempphizandzeta );
    k2scal[1] = tempphizandzeta[0];
    k2scal[2] = tempphizandzeta[1];
    
    
    //calculate k3
    ScalarMultiply( deltat*0.5 , k2vec[0] , temps1 );
    SumVectors( temps1 , s1 , temps1 );
    ScalarMultiply( deltat*0.5 , k2vec[1] , temps2 );
    SumVectors( temps2 , s2 , temps2 );
    ScalarMultiply( deltat*0.5 , k2vec[2] , tempLhat );
    SumVectors( tempLhat , Lhat , tempLhat );
    tempy = y + 0.5*deltat*k2scal[0];
    CalcAngularMomentaDerivs( M , kappa1 , kappa2 , m1 , m2 , tempy , temps1 , temps2 , tempLhat , k3vec[0] , k3vec[1] , k3vec[2] );
    k3scal[0] = Calcdydt( tempy , m1*m2 )/M;
    CalcdPhizdtAnddZetadt( M , m1 , m2 , m1*m2 , tempy , temps1 , temps2 , tempLhat , k3vec[2] , tempphizandzeta );
    k3scal[1] = tempphizandzeta[0];
    k3scal[2] = tempphizandzeta[1];
    
    //calculate k4
    ScalarMultiply( deltat , k3vec[0] , temps1 );
    SumVectors( temps1 , s1 , temps1 );
    ScalarMultiply( deltat , k3vec[1] , temps2 );
    SumVectors( temps2 , s2 , temps2 );
    ScalarMultiply( deltat , k3vec[2] , tempLhat );
    SumVectors( tempLhat , Lhat , tempLhat );
    tempy = y + deltat*k3scal[0];
    CalcAngularMomentaDerivs( M , kappa1 , kappa2 , m1 , m2 , tempy , temps1 , temps2 , tempLhat , k4vec[0] , k4vec[1] , k4vec[2] );
    k4scal[0] = Calcdydt( tempy , m1*m2 )/M;
    CalcdPhizdtAnddZetadt( M , m1 , m2 , m1*m2 , tempy , temps1 , temps2 , tempLhat , k4vec[2] , tempphizandzeta );
    k4scal[1] = tempphizandzeta[0];
    k4scal[2] = tempphizandzeta[1];
    
    //use these to calculate the returns
    angularDerivs[0][0] = (deltat/6.0)*( k1vec[0][0] + 2.0*k2vec[0][0] + 2.0*k3vec[0][0] + k4vec[0][0] );
    angularDerivs[0][1] = (deltat/6.0)*( k1vec[0][1] + 2.0*k2vec[0][1] + 2.0*k3vec[0][1] + k4vec[0][1] );
    angularDerivs[0][2] = (deltat/6.0)*( k1vec[0][2] + 2.0*k2vec[0][2] + 2.0*k3vec[0][2] + k4vec[0][2] );
    angularDerivs[1][0] = (deltat/6.0)*( k1vec[1][0] + 2.0*k2vec[1][0] + 2.0*k3vec[1][0] + k4vec[1][0] );
    angularDerivs[1][1] = (deltat/6.0)*( k1vec[1][1] + 2.0*k2vec[1][1] + 2.0*k3vec[1][1] + k4vec[1][1] );
    angularDerivs[1][2] = (deltat/6.0)*( k1vec[1][2] + 2.0*k2vec[1][2] + 2.0*k3vec[1][2] + k4vec[1][2] );
    angularDerivs[2][0] = (deltat/6.0)*( k1vec[2][0] + 2.0*k2vec[2][0] + 2.0*k3vec[2][0] + k4vec[2][0] );
    angularDerivs[2][1] = (deltat/6.0)*( k1vec[2][1] + 2.0*k2vec[2][1] + 2.0*k3vec[2][1] + k4vec[2][1] );
    angularDerivs[2][2] = (deltat/6.0)*( k1vec[2][2] + 2.0*k2vec[2][2] + 2.0*k3vec[2][2] + k4vec[2][2] );
    scalarDerivs[0] = (deltat/6.0)*( k1scal[0] + 2.0*k2scal[0] + 2.0*k3scal[0] + k4scal[0] );
    scalarDerivs[1] = (deltat/6.0)*( k1scal[1] + 2.0*k2scal[1] + 2.0*k3scal[1] + k4scal[1] );
    scalarDerivs[2] = (deltat/6.0)*( k1scal[2] + 2.0*k2scal[2] + 2.0*k3scal[2] + k4scal[2] );
}


long double NUMCalcThetaL( long double Lhat[3] , long double s1[3] , long double s2[3] , long double m1 , long double m2 , long double eta , long double y )
{
    long double J[3];
    J[0] = eta*Lhat[0]/y + m1*s1[0] + m2*s2[0];
    J[1] = eta*Lhat[1]/y + m1*s1[1] + m2*s2[1];
    J[2] = eta*Lhat[2]/y + m1*s1[2] + m2*s2[2];
    Normalize(J);
    return acosl( DotProduct( J , Lhat ) );
}



//Function to numerically evolve the system
//Inputs:
//      -the initial orbital frequency of the system
//      -the (approximate) final orbital frequency of the system
//      -the mass ratio
//      -the mass normalized spin vector of the first body
//      -the mass normalized spin vector of the second body
//      -the angle between L and s1
//      -the angle between L and s2
//      -the angle between s1 and s2
//      -the quadrupole parameter of body 1
//      -the quadrupole parameter of body 2
//      -the number of points to use in the time evolution
struct Waveform* NumericalEvolver( long double freqInit , long double freqFinal , long double M , long double q , long double chi1 , long double chi2 , long double theta1 , long double theta2 , long double theta12 , long double kappa1 , long double kappa2 , int points )
{
    //initialize the constants of the system
    long double m1 = 1.0 / ( 1.0 + q );
    long double m2 = q*m1;
    long double deltamu = m1 - m2;
    long double eta = m1*m2;
    long double S1 = chi1*m1*m1;
    long double S2 = chi2*m2*m2;
    long double yInit = powl( M*M_PI*freqInit/(300000000) , 1.0/3.0 );
    long double yFinal = powl( M*M_PI*freqFinal/(300000000) , 1.0/3.0 );
    long double tFinal = (M*5.0/(eta*256.0))*( powl( yInit , -8.0 ) -  powl( yFinal , -8.0 ) );
    
    //initialize the dynamical variables
    long double *tValues , *y , *phiz , *zeta , *thetaL;
    tValues = malloc( points*sizeof(long double) );
    phiz = malloc( points*sizeof(long double) );
    zeta = malloc( points*sizeof(long double) );
    thetaL = malloc( points*sizeof(long double) );
    y = malloc( points*sizeof(long double) );
    long double **s1 = (long double **)malloc(points * sizeof(long double*));
    for(int i = 0; i < points; i++)
        s1[i] = (long double *)malloc(3 * sizeof(long double));
    long double **s2 = (long double **)malloc(points * sizeof(long double*));
    for(int i = 0; i < points; i++)
        s2[i] = (long double *)malloc(3 * sizeof(long double));
    long double **Lhat = (long double **)malloc(points * sizeof(long double*));
    for(int i = 0; i < points; i++)
        Lhat[i] = (long double *)malloc(3 * sizeof(long double));
    
    for( int i = 0; i < points; i++ )
    {
        tValues[i] = i*tFinal/(points-1);
    }
    y[0] = yInit;
    phiz[0] = 0.0;
    zeta[0] = 0.0;
    Lhat[0][0] = 0.0;
    Lhat[0][1] = 0.0;
    Lhat[0][2] = 1.0;
    s1[0][0] = S1*sinl( theta1 )/m1;
    s1[0][1] = 0.0;
    s1[0][2] = S1*cosl( theta1 )/m1;
    s2[0][0] = (S2/m2)*( sinl(theta2)*cosl(theta12) );
    s2[0][1] = (S2/m2)*( sinl(theta2)*sinl(theta12) );
    s2[0][2] = (S2/m2)*( cosl(theta2) );
    thetaL[0] = NUMCalcThetaL( Lhat[0] , s1[0] , s2[0] , m1 , m2 , eta , y[0] );


    printf("chi1: %Lf,%Lf,%Lf\t chi2: %Lf,%Lf,%Lf\t, m1:%Lf,m2:%Lf,y0: %Lf\t thetaL:%Lf\n", s1[0][0]/m1,s1[0][1],s1[0][2]/m1, s2[0][0]/m2,s2[0][1]/m2,s2[0][2]/m2, m1, m2, y[0], thetaL[0]);

    printf("Input--> s1: %Lf,%Lf,%Lf\t s2: %Lf,%Lf,%Lf\t, m1:%Lf,m2:%Lf,y0: %Lf\t thetaL:%Lf\n", s1[0][0],s1[0][1],s1[0][2], s2[0][0],s2[0][1],s2[0][2], m1, m2, y[0], thetaL[0]);


    //actually evolve the system using RK4
    long double temp1[3][3] , temp2[3];
    for( int i = 1; i < points; i++ )
    {
        RK4evolver( M , kappa1 , kappa2 , m1 , m2 , y[i-1] , s1[i-1] , s2[i-1] , Lhat[i-1] , tValues[i]-tValues[i-1] , temp1 , temp2 );
        y[i] = y[i-1] + temp2[0];
        phiz[i] = phiz[i-1] + temp2[1];
        zeta[i] = zeta[i-1] + temp2[2];
        SumVectors( s1[i-1] , temp1[0] , s1[i] );
        SumVectors( s2[i-1] , temp1[1] , s2[i] );
        SumVectors( Lhat[i-1] , temp1[2] , Lhat[i] );
        
        thetaL[i] = NUMCalcThetaL( Lhat[i] , s1[i] , s2[i] , m1 , m2 , eta , y[i] );
    }
    
    
    
    //output the relevant data in csv form to a text document for plotting purposes
    
    char filename[100];
    sprintf(filename,"plotNumerical - k1 = %.2Lf , k2 = %.2Lf.txt" , kappa1 , kappa2);
    FILE *filey = fopen(filename, "w");
    fprintf(filey, "value number,time,chieff,deltachi,J,phiz,zeta,thetaL \n");
    long double deltachi , chieff , J;
    for( int i = 0; i < points; i++ )
    {
        deltachi = CalcDeltachi( s1[i] , s2[i] , Lhat[i] );
        chieff = CalcChieff( s1[i] , s2[i] , Lhat[i] );
        J = CalcJ( s1[i] , s2[i] , Lhat[i] , m1 , m2 , y[i] , eta );
        fprintf(filey, "val#%i,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le \n", i ,  tValues[i], chieff , deltachi , J , phiz[i] , zeta[i] , thetaL[i], DotProduct( s1[i] , Lhat[i] ), DotProduct( s2[i] , Lhat[i] ) );
   }
    fflush(filey);
    
    
    
    struct Waveform* returner = GenWaveform( M , tValues , y , phiz , zeta , thetaL , points , 0 , 0 , 0 , 0 , eta );
    free(tValues);
    free(phiz);
    free(zeta);
    free(thetaL);
    free(y);
    for (int i = 0; i < points; i++)
        free(s1[i]);
    free(s1);
    for (int i = 0; i < points; i++)
        free(s2[i]);
    free(s2);
    for (int i = 0; i < points; i++)
        free(Lhat[i]);
    free(Lhat);
    return returner;
    
}
