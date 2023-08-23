//
//  BHNSAnalyticEvolver.c
//  
//
//  Created by Michael Lahaye on 2022-04-02.
//

#include "BHNSAnalyticalEvolver.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "BHNSVectors.h"
#include "BHNSNumericalEvolver.h"
#include "BHNSWaveform.h"
#include <complex.h>



//function to calculate the sign of a number, used to find the sign of the initial derivatives to calculate the relative sign between the amplitudes
//Inputs:
//      -the number
//Outputs:
//      -the sign
long double sign( long double x )
{

    if(x<0)
    {
        return -1.0;
    }
    if(x==0.0)
    {
        return 0.0;
    }
    return 1.0;

}

void sortRoots( long double roots[3] )
{
    long double temp;
    if( fabsl(roots[0]) > fabsl(roots[1]) )
    {
        temp = roots[0];
        roots[0] = roots[1];
        roots[1] = temp;
    }
    if( fabsl(roots[1]) > fabsl(roots[2]) )
    {
        temp = roots[1];
        roots[1] = roots[2];
        roots[2] = temp;
    }
    if( fabsl(roots[0]) > fabsl(roots[1]) )
    {
        temp = roots[0];
        roots[0] = roots[1];
        roots[1] = temp;
    }
}



//returns the plus, minus, and third root in that order.
void cubicRootFinder( long double x3coeff , long double x2coeff , long double x1coeff , long double x0coeff , long double returnVec[3] )
{

    long double a = x2coeff/x3coeff;
    long double b = x1coeff/x3coeff;
    long double c = x0coeff/x3coeff;
    
    long double Q = ( a*a - 3.0*b )/9.0;
    long double R = ( 2.0*a*a*a - 9.0*a*b + 27.0*c )/54.0;
    long double S = Q*Q*Q - R*R;
    
    if( S>0 )
    {
        long double phi = (1.0/3.0)*acosl( R/powl( Q , 1.5 ) );
        returnVec[0] = -2.0*powl( Q , 0.5 )*cosl( phi - 2.0*M_PI/3.0 ) - a/3.0;
        returnVec[1] = -2.0*powl( Q , 0.5 )*cosl( phi ) - a/3.0;
        returnVec[2] = -2.0*powl( Q , 0.5 )*cosl( phi + 2.0*M_PI/3.0 ) - a/3.0;
        sortRoots( returnVec );
        return;
    }
    if( S==0 )
    {
        returnVec[0] = sign(R)*powl( Q , 0.5 ) - a/3;
        returnVec[1] = -2.0*sign(R)*powl( Q , 0.5 ) - a/3.0;
        returnVec[2] = sign(R)*powl( Q , 0.5 ) - a/3;
        sortRoots( returnVec );
        return;
    }
    
    returnVec[0] = 0;
    returnVec[1] = 0;
    returnVec[2] = 0;
}



//Function to calculate the initial sign of the derivative of deltachi
//Inputs:
//      -the quadrupole parameter of body 1
//      -the quadrupole parameter of body 2
//      -the mass of body 1
//      -the mass of body 2
//      -the initial value of the PN parameter y
//      -the initial mass normalized spin vector s1
//      -the initial mass normalized spin vector s2
//      -the initial normalized orbital angular momentum vector Lhat
//Outputs:
//      -the sign
long double CalcdeltachiAMPSign( long double M , long double kappa1 , long double kappa2 , long double m1 , long double m2 , long double y , long double s1[3] , long double s2[3] , long double Lhat[3] )
{

    long double ds1[3] , ds2[3] , dLhat[3];
    CalcAngularMomentaDerivs( M , kappa1 , kappa2 , m1 , m2 , y , s1 , s2 , Lhat , ds1 , ds2 , dLhat );
    long double ddeltachi = DotProduct( ds1 , Lhat ) - DotProduct( ds2 , Lhat ) + DotProduct( s1 , dLhat ) - DotProduct( s2 , dLhat );
    
    return sign( ddeltachi );
    
}



//Function to calculate the initial sign of the derivative of chieff
//Inputs:
//      -the quadrupole parameter of body 1
//      -the quadrupole parameter of body 2
//      -the mass of body 1
//      -the mass of body 2
//      -the initial value of the PN parameter y
//      -the initial mass normalized spin vector s1
//      -the initial mass normalized spin vector s2
//      -the initial normalized orbital angular momentum vector Lhat
//Outputs:
//      -the sign
long double CalcchieffAMPSign( long double M , long double kappa1 , long double kappa2 , long double m1 , long double m2 , long double y , long double s1[3] , long double s2[3] , long double Lhat[3] )
{

    long double ds1[3] , ds2[3] , dLhat[3];
    CalcAngularMomentaDerivs( M , kappa1 , kappa2 , m1 , m2 , y , s1 , s2 , Lhat , ds1 , ds2 , dLhat );
    long double dchieff = DotProduct( ds1 , Lhat ) + DotProduct( ds2 , Lhat ) + DotProduct( s1 , dLhat ) + DotProduct( s2 , dLhat );
    return sign( dchieff );
    
}



//Function to calculate the coefficients B0 and B1
//Inputs:
//      -the total angular momentum J
//      -the orbital angular momentum L
//      -the magnitude of body 1's spin
//      -the magnitude of body 2's spin
//      -the current value of the PN parameter
//      -the symmetric mass ratio
//      -the mass difference
//      -the vector to return the values of B
//Outputs:
//      -none
void CalcB( long double J , long double L , long double S1 , long double S2 , long double y , long double eta , long double deltamu , long double returnvector[2] )
{
    returnvector[0] = (0.5*y/(eta*eta))*( -2.0*eta*( J*J - L*L ) + deltamu*( S1*S1 - S2*S2 ) - deltamu*deltamu*( 2*L*L + S1*S1 + S2*S2 ) );
    returnvector[1] = 1.0;
}



//Function to calculate the coefficients C0, C1, and C2
//Inputs:
//      -the total angular momentum J
//      -the orbital angular momentum L
//      -the magnitude of body 1's spin
//      -the magnitude of body 2's spin
//      -the current value of the PN parameter
//      -the symmetric mass ratio
//      -the mass difference
//      -the vector to return the values of B
//Outputs:
//      -none
void CalcC( long double J , long double L , long double S1 , long double S2 , long double y , long double eta , long double deltamu , long double returnvector[3] )
{
    returnvector[0] = 2.0*deltamu*( J*J - L*L - S1*S1 - S2*S2 )/eta;
    returnvector[1] = 0.5*y*( ( 1 + deltamu*deltamu )*( S1*S1 - S2*S2 ) - 2.0*deltamu*( 2.0*L*L + S1*S1 + S2*S2 ) )/(eta*eta);
    returnvector[2] = -deltamu;
}



//Function to calculate the coefficients D0, D1, D2, and D3
//Inputs:
//      -the total angular momentum J
//      -the orbital angular momentum L
//      -the magnitude of body 1's spin
//      -the magnitude of body 2's spin
//      -the current value of the PN parameter
//      -the symmetric mass ratio
//      -the mass difference
//      -the vector to return the values of B
//Outputs:
//      -none
void CalcD( long double J , long double L , long double S1 , long double S2 , long double y , long double eta , long double deltamu , long double returnvector[4] )
{
    returnvector[0] = -y*( J*J - L*L - S1*S1 - 2.0*S1*S2 - S2*S2 )*( J*J - L*L - S1*S1 + 2.0*S1*S2 - S2*S2 )/(eta*eta);
    returnvector[1] = 2.0*( J*J - L*L - S1*S1 - S2*S2 )/eta;
    returnvector[2] = 0.5*y*( 2.0*eta*( J*J - L*L ) - ( 2*L*L + S1*S1 + S2*S2 ) + deltamu*( S1*S1 - S2*S2 ) )/(eta*eta);
    returnvector[3] = -1.0;
}



//Function to compute the first three derivatives of deltachi at the initial time.  For convenience of writing, we'll relabel deltachi as X, relabel chieff and Y, and (dX/dt)^2 = (9/4)y^11*Q.  This is just to simplify the code.
//Inputs:
//      -the A coefficients for deltachi
//      -the A coefficients for chieff
//      -the cubic coefficients B
//      -the cubic coefficients C
//      -the cubic coefficients D
//      -deltachi
//      -chieff
//      -the current value of the PN parameter
//      -the mass difference
//      -the sign of the initial derivative of deltachi
//      -the sign of the initial derivative of chieff
//      -the location to return the first three derivatives
//Outputs:
//      -none
void CalcInitDeltachiDerivs( long double AX[2] , long double AY[2] , long double B[2] , long double C[3] , long double D[4] , long double X , long double Y , long double y , long double deltamu , long double Xsign , long double Ysign , long double returnVec[3] )
{
    
    long double AXsum = 1 + y*AX[0]*X + y*AX[1]*Y;
    long double AYsum = AY[0]*X + AY[1]*Y;
    long double Q = AXsum*AXsum*( deltamu*X*X*X + (B[0]+B[1]*Y)*X*X + (C[0]+C[1]*Y+C[2]*Y*Y)*X + (D[0]+D[1]*Y+D[2]*Y*Y+D[3]*Y*Y*Y) );
    
    long double dQdX = 2.0*y*AX[0]*Q/AXsum + ( 1 + y*AX[0]*X + y*AX[1]*Y )*( 1 + y*AX[0]*X + y*AX[1]*Y )*( 3.0*deltamu*X*X + 2.0*(B[0]+B[1]*Y)*X + (C[0]+C[1]*Y+C[2]*Y*Y) ) ;
    long double dQdY = 2.0*y*AX[1]*Q/AXsum + ( 1 + y*AX[0]*X + y*AX[1]*Y )*( 1 + y*AX[0]*X + y*AX[1]*Y )*( B[1]*X*X + (C[1]+2.0*C[2]*Y)*X + (D[1]+2.0*D[2]*Y+3.0*D[3]*Y*Y) );
    long double d2QdX2 = 2.0*y*y*AX[0]*AX[0]*Q/( AXsum*AXsum ) + 4.0*AXsum*y*AX[0]*( 3.0*deltamu*X*X + 2.0*(B[0]+B[1]*Y)*X + (C[0]+C[1]*Y+C[2]*Y*Y) ) + AXsum*AXsum*( 2.0*B[1]*Y + 2.0*B[0] + 6.0*deltamu*X );
    long double d2QdXdY = 2*y*y*AX[0]*AX[1]*Q/( AXsum*AXsum ) + 2.0*AXsum*y*AX[1]*( 3.0*deltamu*X*X + 2.0*(B[0]+B[1]*Y)*X + (C[0]+C[1]*Y+C[2]*Y*Y) ) + 2.0*AXsum*y*AX[0]*( B[1]*X*X + (2.0*C[2]*Y+C[1])*X + D[1] + 2.0*D[2]*Y + 3.0*D[3]*Y*Y ) + AXsum*AXsum*( 2.0*B[1]*X + 2*C[2]*Y + C[1] );
    long double d2QdY2 = 2.0*y*y*AX[1]*AX[1]*Q/( AXsum*AXsum ) + 4.0*AXsum*y*AX[1]*( B[1]*X*X + (C[1]+2.0*C[2]*Y)*X + (D[1]+2.0*D[2]*Y+3.0*D[3]*Y*Y) ) + AXsum*AXsum*( 2.0*C[2]*X + 6.0*D[3]*Y + 2.0*D[2] );
    
    long double ratio = fabsl( y*AYsum/AXsum );
    long double dratiodX = fabsl( y*AY[0]/AXsum - AYsum*y*AX[0]/( AXsum*AXsum ) );
    long double dratiodY = fabsl( y*AY[1]/AXsum - AYsum*y*AX[1]/( AXsum*AXsum ) );
    
    long double dXdt = Xsign*3.0*powl( y , 5.5 )*powl( Q , 0.5 )/2.0;
    long double dYdt = (Ysign/Xsign)*ratio*dXdt;
    returnVec[0] = dXdt;
    returnVec[1] = (9.0/4.0)*powl( y , 11.0 )*0.5*( dQdX + Xsign*Ysign*ratio*dQdY );
    returnVec[2] = (9.0/8.0)*powl( y , 11.0 )*( dXdt*d2QdX2 + 2.0*dYdt*d2QdXdY + dYdt*dYdt*d2QdY2/dXdt + Xsign*Ysign*dQdY*( dXdt*dratiodX + dYdt*dratiodY ) );
    
}



//using the initial values of chieff and deltachi, compute the first three derivatives of deltachi at the initial time.  From these, compute the initial average of deltachi and the initial phase of deltachi in its precession.  Then use these and the initial value of chieff to calculate the average of chieff.
//Inputs:
//      -the A coefficients for deltachi
//      -the A coefficients for chieff
//      -the cubic coefficients B
//      -the cubic coefficients C
//      -the cubic coefficients D
//      -deltachi
//      -chieff
//      -the current value of the PN parameter
//      -the mass difference
//      -the sign of the initial derivative of deltachi
//      -the sign of the initial derivative of chieff
//      -the location to return the first three derivatives
//Outputs:
//      -none
void CalcInitAVGandpsi( long double Adeltachi[2] , long double Achieff[2] , long double B[2] , long double C[3] , long double D[4] , long double deltachi , long double chieff , long double y , long double deltamu , long double deltachisign , long double chieffsign , long double deltachiAVG[1] , long double chieffAVG[1] , long double psi[1] )
{
    
    long double temp[3] , AMP , dphidt;
    CalcInitDeltachiDerivs( Adeltachi , Achieff , B , C , D , deltachi , chieff , y , deltamu , deltachisign , chieffsign , temp );
    
    //using the approximate sin relation for deltachi, we have deltachiAVG(0) = deltachi(0) - secondderiv(0)*firstderiv(0)/thirdderiv(0)
    deltachiAVG[0] = deltachi - temp[1]*temp[0]/temp[2];
    
    //using the approximate sin relation for deltachi, the initial phase is tan(psi(0)) = -| ( secondderiv(0)/thirdderiv(0) )*√( thirdderiv(0)/firstderiv(0) ) |
    long double tan = fabsl( (temp[1]/temp[2])*powl( fabsl(temp[2]/temp[0]) , 0.5 ) );
    //we then have to use the signs of the first and second derivative to figure out which quadrant the angle is currently in and adjust accordingly
    psi[0] = atanl(tan);
    if( sign(temp[1])>0 && sign(temp[0])>0 )
        psi[0] = -psi[0];
    if( sign(temp[1])<0 && sign(temp[0])<0 )
        psi[0] = M_PI - psi[0];
    if( sign(temp[1])>0 && sign(temp[0])<0 )
        psi[0] = M_PI + psi[0];
    dphidt = powl( -temp[2]/temp[0] , 0.5 );
    
    //AMP = ( deltachi - deltachiAVG[0])/sinl(psi[0]);
    AMP = powl( (temp[0]*temp[0])/(dphidt*dphidt) + (temp[1]*temp[1])/(dphidt*dphidt*dphidt*dphidt) , 0.5);
    
    //chieffAVG[0] = chieff - y*( Achieff[0]*deltachiAVG[0] + Achieff[1]*chieff )*AMP*( sinl(psi[0])*sinl(psi[0]) )/( 1 - y*chieff ) - 0.5*y*Achieff[0]*AMP*AMP*( sinl(psi[0])*sinl(psi[0])*sinl(psi[0])*sinl(psi[0]) )/( 1 - y*chieff );
    
    /*
    long double m, rooties[3];
    cubicRootFinder( deltamu , B[0]+chieff*B[1] , C[0]+chieff*C[1]+chieff*chieff*C[2] , D[0]+chieff*D[1]+chieff*chieff*D[2]+chieff*chieff*chieff*D[3] , rooties );
    m = (rooties[0]-rooties[1])/(rooties[2]-rooties[1]);
    */
    long double psibar = 0.5*psi[0] + 0.25*M_PI;
    chieffAVG[0] = chieff - 2.0*y*( Achieff[0]*(deltachiAVG[0]-AMP) + Achieff[1]*chieff )*AMP*( sinl(psibar)*sinl(psibar) - 0.5 )/( 1 - y*chieff ) - 2.0*y*Achieff[0]*AMP*AMP*( sinl(psibar)*sinl(psibar)*sinl(psibar)*sinl(psibar) - 0.5 )/( 1 - y*chieff );
    
}



//calculate the amplitudes of deltachi and chieff, and the derivative of the phase psi
void CalcAMPsAndPhaseDeriv( long double Adeltachi[2] , long double Achieff[2] , long double B[2] , long double C[3] , long double D[4] , long double deltamu , long double y , long double deltachiAVG , long double chieffAVG , long double deltachisign , long double chieffsign , long double returnVec[3] )
{
    long double x3coeff , x2coeff , x1coeff , x0coeff , N0 , N1 , ratio , roots[3] , AAVG , AAMP;
    N1 = (deltachisign*chieffsign)*fabsl( ( y*Achieff[0]*deltachiAVG + y*Achieff[1]*chieffAVG )/( 1.0 + y*Adeltachi[0]*deltachiAVG + y*Adeltachi[1]*chieffAVG ) );
    N0 = chieffAVG - N1*deltachiAVG;
    
    x0coeff = D[3]*N0*N0*N0 + D[2]*N0*N0 + D[1]*N0 + D[0];
    x1coeff = 3.0*D[3]*N0*N0*N1 + C[2]*N0*N0 + 2.0*D[2]*N0*N1 + C[1]*N0 + D[1]*N1 + C[0];
    x2coeff = 3.0*D[3]*N0*N1*N1 + 2.0*C[2]*N0*N1 + D[2]*N1*N1 + B[1]*N0 + C[1]*N1 + B[0];
    x3coeff = D[3]*N1*N1*N1 + C[2]*N1*N1 + B[1]*N1 + deltamu;
    
    cubicRootFinder( x3coeff , x2coeff , x1coeff , x0coeff , roots );
    returnVec[0] = 0.5*fabsl( roots[1] - roots[0] );
    returnVec[1] = N1*returnVec[0];
    
    AAVG = 1.0 + y*Adeltachi[0]*deltachiAVG + y*Adeltachi[1]*chieffAVG;
    AAMP = y*Adeltachi[0]*returnVec[0] + y*Adeltachi[1]*returnVec[1];
    returnVec[2] = powl( (9.0/4.0)*x3coeff*powl(y,11.0)*( 0.5*( roots[2] - deltachiAVG )*( 2*AAVG*AAVG + AAMP*AAMP ) - returnVec[0]*AAMP*AAVG ) , 0.5 );
}


//returns the derivative of JAVG and then JAMPSIN and JAMPCOS in that order to return vector.
void CalcdJAVGAndJAMPs( long double dydt , long double dpsidt , long double L , long double y , long double deltamu , long double JAVG , long double deltachiAVG , long double deltachiAMP , long double chieffAVG , long double chieffAMP , long double returnVec[3] )
{
    
    long double QB , QC;
    QB = (0.5*L/(y*JAVG*JAVG))*dydt*( 2*L + chieffAVG + deltamu*deltachiAVG );
    QC = ( 0.5*L/(JAVG*y) )*dydt*( chieffAMP + deltamu*deltachiAMP );
    
    returnVec[1] = QB*QC/( dpsidt*dpsidt + QB*QB );
    returnVec[2] = dpsidt*returnVec[1]/QB;
    returnVec[0] = (-QB*JAVG + 0.5*QC*returnVec[1]/(JAVG));
    
}



void CalcDeltachiAVGDeriv( long double B[2] , long double C[3] , long double D[4] , long double Adeltachi[2] , long double Achieff[2] , long double deltamu , long double eta , long double S1 , long double S2 , long double y , long double dydt , long double JAVG , long double dJAVGdt , long double L , long double dLdt , long double deltachiplus , long double deltachiminus , long double chieffplus , long double chieffminus , long double tempInit[2] )
{
    long double dB0dt , dC0dt , dC1dt , dD0dt , dD1dt , dD2dt, plusratio , minusratio , LHSplus , LHSminus , coeff1plus , coeff2plus , coeff1minus , coeff2minus ;
    
    dB0dt = dydt*B[0]/y - 2.0*y*( eta*( JAVG*dJAVGdt - L*dLdt ) + deltamu*deltamu*L*dLdt )/(eta*eta);
    dC0dt = 4.0*deltamu*( JAVG*dJAVGdt - L*dLdt )/eta;
    dC1dt = dydt*C[1]/y - 4.0*y*deltamu*L*dLdt/(eta*eta);
    dD0dt = dydt*D[0]/y - 4.0*y*( JAVG*dJAVGdt - L*dLdt )*( JAVG*JAVG - L*L - S1*S1 - S2*S2 )/(eta*eta);
    dD1dt = 4.0*( JAVG*dJAVGdt - L*dLdt )/eta;
    dD2dt = dydt*D[2]/y + 2.0*y*( eta*( JAVG*dJAVGdt - L*dLdt ) - L*dLdt )/(eta*eta);
    
    plusratio = ( Achieff[0]*y*deltachiplus + Achieff[1]*y*chieffplus )/( 1 + Adeltachi[0]*y*deltachiplus + Adeltachi[1]*y*chieffplus);
    minusratio = ( Achieff[0]*y*deltachiminus + Achieff[1]*y*chieffminus )/( 1 + Adeltachi[0]*y*deltachiminus + Adeltachi[1]*y*chieffminus);
    
    LHSplus = -( dB0dt*deltachiplus*deltachiplus + ( dC0dt + dC1dt*chieffplus )*deltachiplus + ( dD0dt + dD1dt*chieffplus + dD2dt*chieffplus*chieffplus ) );
    LHSminus = -( dB0dt*deltachiminus*deltachiminus + ( dC0dt + dC1dt*chieffminus )*deltachiminus + ( dD0dt + dD1dt*chieffminus + dD2dt*chieffminus*chieffminus ) );
    
    coeff1plus = 3.0*deltamu*deltachiplus*deltachiplus + 2.0*( B[0] + B[1]*chieffplus )*deltachiplus + C[0] + C[1]*chieffplus + C[2]*chieffplus*chieffplus;
    coeff1minus = 3.0*deltamu*deltachiminus*deltachiminus + 2.0*( B[0] + B[1]*chieffminus )*deltachiminus + C[0] + C[1]*chieffminus + C[2]*chieffminus*chieffminus;
    
    coeff2plus = B[1]*deltachiplus*deltachiplus + (C[1]+2.0*C[2]*chieffplus)*deltachiplus + D[1] + 2.0*D[2]*chieffplus + 3.0*D[3]*chieffplus*chieffplus;
    coeff2minus = B[1]*deltachiminus*deltachiminus + (C[1]+2.0*C[2]*chieffminus)*deltachiminus + D[1] + 2.0*D[2]*chieffminus + 3.0*D[3]*chieffminus*chieffminus;
    
    long double dXp , dXm , dYp , dYm;
    
    dXp = LHSplus/( coeff1plus + plusratio*coeff2plus );
    dYp = dXp*plusratio;
    dXm = LHSminus/( coeff1minus + minusratio*coeff2minus );
    dYm = dXm*minusratio;
    
    tempInit[0] = 0.5*( dXp + dXm );
    tempInit[1] = 0.5*( dYp + dYm );
}


void CalcInts( long double a , long double b , long double Int[5] )
{
    long double aprime , bprime , commondenominator;
    aprime = powl( 1.0 - a*a , 0.5 );
    bprime = powl( 1.0 - b*b , 0.5 );
    commondenominator = 1.0/( aprime*bprime*(a-b) );

    Int[0] = ( bprime*a - aprime*b )*commondenominator;
    Int[1] = ( aprime - bprime )*commondenominator;
    Int[2] = ( ( (a-b)*bprime - a )*aprime + b*bprime )*commondenominator/( a*b );
    Int[3] = -( ( (a*a-b*b)*bprime - a*a )*aprime + b*b*bprime )*commondenominator/( a*a*b*b );
    Int[4] = 0.5*( ( ( (b*b+2.0)*a*a*a - a*a*b*b*b - 2.0*b*b*b )*bprime - 2*a*a*a )*aprime + 2*bprime*b*b*b )*commondenominator/( a*a*a*b*b*b );
    
}



void CalcdphizAVGAndphizAMP( long double q , long double eta , long double deltamu , long double S1 , long double S2 , long double y , long double JAVG , long double JAMPSIN , long double JAMPCOS , long double deltachiAVG , long double deltachiAMP , long double chieffAVG , long double chieffAMP , long double Achieff[2] , long double kappa1 , long double kappa2 , long double mu1 , long double mu2 , long double dpsi , long double psi , long double returnVector[4] )
{
    long double Q1 , Q2 , Q3 , Q4 , Q5 , Q6 , Q7 , Q8 , Q9 , Q10, Q11 , Q12 , Q13 , Q14 , Q15 , Q16 , Q17 , Del , H0 , H1 , H2 , H3 , Hminus , Hplus , Int[5];
    
    Q1 = 1.0 +1.5*( 1.0 - y*chieffAVG )/eta;
    Q2 = -1.5*y*chieffAMP/eta;
    Q3 = ( 1.0 - Q1 )/( 1 + q );
    Q4 = -Q2/( 1 + q );
    Q5 = 4.0*(1.0-q)*( S1*S1 - S2*S2 ) - (1.0+q)*( deltamu*deltachiAVG + chieffAVG )*( deltamu*deltachiAVG + (1.0-4.0*eta)*chieffAVG );
    Q6 = -(1.0+q)*( ( deltamu*deltachiAVG + chieffAVG )*( deltamu*deltachiAMP + (1.0-4.0*eta)*chieffAMP ) + ( deltamu*deltachiAMP + chieffAMP )*( deltamu*deltachiAVG + (1.0-4.0*eta)*chieffAVG ) );
    Q7 = -(1.0+q)*( deltamu*deltachiAMP + chieffAMP )*( deltamu*deltachiAMP + (1.0-4.0*eta)*chieffAMP );
    Q8 = 2.0*JAVG - deltamu*deltachiAVG - chieffAVG - 2*eta/y;
    Q9 = 4.0*JAVG - Q8;
    Q10 = deltamu*deltachiAMP + chieffAMP;
    
    Q11 = -3.0*y*( JAVG - eta/y )*( JAVG + eta/y )*( (kappa2-1)*mu1*( chieffAVG - deltachiAVG ) + (kappa1-1)*mu2*( chieffAVG + deltachiAVG ) )/eta + 3.0*y*( S1 - S2 )*( S1 + S2 )*( (kappa2-1)*mu1*( chieffAVG - deltachiAVG ) - (kappa1-1)*mu2*( chieffAVG + deltachiAVG ) )/eta;
    Q12 = -3.0*y*( JAVG - eta/y )*( JAVG + eta/y )*( (kappa2-1)*mu1*( chieffAMP - deltachiAMP ) + (kappa1-1)*mu2*( chieffAMP + deltachiAMP ) )/eta + 3.0*y*( S1 - S2 )*( S1 + S2 )*( (kappa2-1)*mu1*( chieffAMP - deltachiAMP ) - (kappa1-1)*mu2*( chieffAMP + deltachiAMP ) )/eta;
    Q13 = 1.5*( chieffAVG + deltamu*deltachiAVG );
    Q14 = 1.5*( chieffAMP + deltamu*deltachiAMP );
    Q15 = 4.0*y*Achieff[0]*deltachiAVG*deltachiAVG + 2.0*( 4.0*y*Achieff[1]*chieffAVG + mu2*kappa1 - mu1*kappa2 + deltamu )*deltachiAVG + 4.0*y*Achieff[0]*chieffAVG*chieffAVG + 2.0*( mu2*kappa1 + mu1*kappa2 - 1 )*chieffAVG;
    Q16 = 2.0*( 4.0*( chieffAVG*Achieff[1] + Achieff[0]*deltachiAVG )*y - mu1*kappa2 + mu2*kappa1 + deltamu )*deltachiAMP + 2.0*( 4.0*( chieffAVG*Achieff[1] + Achieff[0]*deltachiAVG )*y + mu1*kappa2 + mu2*kappa1 - 1 )*chieffAMP;
    Q17 = 4.0*y*( ( chieffAMP*chieffAMP + deltachiAMP*deltachiAMP )*Achieff[0] + 2.0*Achieff[1]*chieffAMP*deltachiAMP );

    
    
    H0 = ( Q3*Q5 + Q11 + Q13*Q15 )/( Q8*Q9 );
    H1 = ( Q3*Q6 + Q5*Q4 + Q12 + Q13*Q16 + Q14*Q15 )/( Q8*Q9 );
    H2 = ( Q3*Q7 + Q4*Q6 + Q13*Q17 + Q14*Q16 )/( Q8*Q9 );
    H3 = ( Q4*Q7 + Q14*Q17 )/( Q8*Q9 );
    Hminus = -Q10/Q8;
    Hplus = Q10/Q9;
    
    
    long double Phi0 = 0.5*JAVG*powl( y , 6.0 )*( Q1 + (  H2*Hplus*Hminus - H3*Hminus - H3*Hplus )/( Hplus*Hplus*Hminus*Hminus ) );
    long double Phis = 0.5*JAVG*powl( y , 6.0 )*( Q2 + H3/( Hplus*Hminus ) );
    long double Phip = 0.5*JAVG*powl( y , 6.0 )*( H0*Hplus*Hplus*Hplus - H1*Hplus*Hplus + H2*Hplus - H3 )/( ( Hplus - Hminus )*Hplus*Hplus*powl( 1 - Hplus*Hplus , 0.5 ) );
    long double Phim = -0.5*JAVG*powl( y , 6.0 )*( H0*Hminus*Hminus*Hminus - H1*Hminus*Hminus + H2*Hminus - H3 )/( ( Hplus - Hminus )*Hminus*Hminus*powl( 1 - Hminus*Hminus , 0.5 ) );
    long double Theta0 = (2.0*eta/y + deltamu*deltachiAVG + chieffAVG)/( 2.0*JAVG );
    long double Thetas = (deltamu*deltachiAMP + chieffAMP)/( 2.0*JAVG );
    long double psi0 = fmod( psi , 2*M_PI );
    if(psi0 > M_PI)
    {
        psi0 = psi0 - 2*M_PI;
    }
    if(psi0 < -M_PI)
    {
        psi0 = psi0 + 2*M_PI;
    }
    long double tempplus = 2*atanl( ( tanl( psi0/2 ) + Hplus )/powl( 1 - Hplus*Hplus , 0.5 ) ) - psi0 - asinl( Hplus );
    long double tempminus = 2*atanl( ( tanl( psi0/2 ) + Hminus )/powl( 1 - Hminus*Hminus , 0.5 ) ) - psi0 - asinl( Hminus );
    
    returnVector[0] = Phi0 + Phip + Phim;
    returnVector[1] = ( Phip*tempplus + Phim*tempminus - Phis*cosl( psi0 ) )/dpsi;
    returnVector[2] = -Theta0*returnVector[0] - Thetas*Phis/2 - Thetas*Phip*(powl( 1 - Hplus*Hplus ,0.5 )-1)/Hplus  - Thetas*Phim*(powl( 1 - Hminus*Hminus ,0.5 )-1)/Hminus ;
    returnVector[3] = -Theta0*returnVector[1] + ( - Thetas*Phi0*cosl(psi) + 0.25*Thetas*Phis*sinl( 2*psi ) + Thetas*Phip*tempplus/Hplus + Thetas*Phim*tempminus/Hminus  )/dpsi;
    
    /*
    Q1 = 4.0*(1.0-q)*( S1*S1 - S2*S2 ) - (1.0+q)*( deltamu*deltachiAVG + chieffAVG )*( deltamu*deltachiAVG + (1.0-4.0*eta)*chieffAVG );
    Q2 = -(1.0+q)*( ( deltamu*deltachiAVG + chieffAVG )*( deltamu*deltachiAMP + (1.0-4.0*eta)*chieffAMP ) + ( deltamu*deltachiAMP + chieffAMP )*( deltamu*deltachiAVG + (1.0-4.0*eta)*chieffAVG ) );
    Q3 = -(1.0+q)*( deltamu*deltachiAMP + chieffAMP )*( deltamu*deltachiAMP + (1.0-4.0*eta)*chieffAMP );
    Q4 = 1.0 +1.5*( 1.0 - y*chieffAVG )/eta;
    Q5 = -1.5*y*chieffAMP/eta;
    Q6 = -1.5*(1.0+q)*( 1.0 - y*chieffAVG )/q;
    Q7 = 1.5*(1.0+q)*y*chieffAMP/q;
    Q8 = 2.0*JAVG - deltamu*deltachiAVG - chieffAVG - 2*eta/y;
    Q9 = JAMPCOS;
    Q10 = JAMPSIN - deltamu*deltachiAMP - chieffAMP;
    Q11 = 2.0*JAVG + deltamu*deltachiAVG + chieffAVG + 2*eta/y;
    Q13 = JAMPSIN + deltamu*deltachiAMP + chieffAMP;
    
    H0 = Q6*Q1;
    H1 = ( Q2*Q6 + Q1*Q7 );
    H2 = ( Q3*Q6 + Q2*Q7 );
    H3 = Q3*Q7;
    Hminus = Q10/Q8;
    Hplus = Q13/Q11;
    
    returnVector[0] = 0.5*powl( y , 6.0 )*( JAVG*Q4 + H0*JAVG/(Q11*Q8) + 0.5*H2*JAVG/(Q11*Q8) + JAMPSIN*Q5*0.5 + 0.5*H1*JAMPSIN/(Q11*Q8) + 3.0*H3*JAMPSIN/(Q11*Q8*8.0) - 0.5*H1*JAVG*Q10/(Q11*Q8*Q8) - 3.0*H3*JAVG*Q10/(Q11*Q8*Q8*8.0) - H1*JAVG*Q13*0.5/(Q11*Q8*Q11) - 3.0*H3*JAVG*Q13/(Q11*Q8*Q11*8.0) );
    returnVector[1] = JAVG*powl( y , 6.0 )*( 0.5*Q5 + H0*Int[1] + H1*Int[2] + H2*Int[3] + H3*Int[4] );
    */
}



long double ANALCalcThetaL( long double J , long double L , long double chieff , long double deltachi , long double deltamu )
{
    return acosl( 0.5*( 2*L + chieff + deltamu*deltachi )/J );
}



struct Asys* Initsys( int points , long double M , long double freqInit , long double freqFinal , long double q , long double chi1 , long double chi2 , long double theta1 , long double theta2 , long double theta12 , long double kappa1 , long double kappa2 )
{
    struct Asys* returnPTR = malloc( sizeof(struct Asys) );
    
    //initialize all the constants
    returnPTR -> M = M;
    returnPTR -> q = q;
    returnPTR -> m1 = 1.0/( 1.0 + q );
    returnPTR -> m2 = q * ((*returnPTR).m1);
    returnPTR -> dm = (*returnPTR).m1 - (*returnPTR).m2;
    returnPTR -> eta = ((*returnPTR).m1)*((*returnPTR).m2);
    returnPTR -> S1 = chi1*((*returnPTR).m1)*((*returnPTR).m1);
    returnPTR -> S2 = chi2*((*returnPTR).m2)*((*returnPTR).m2);
    returnPTR -> k1 = kappa1;
    returnPTR -> k2 = kappa2;
    returnPTR -> Adc[0] = ( kappa2 - kappa1 )/4.0;
    returnPTR -> Adc[1] = -( kappa1 + kappa2 + 2.0 )/4.0;
    returnPTR -> Ace[0] = ( kappa1 + kappa2 - 2.0 )/4.0;
    returnPTR -> Ace[1] = - ( (*returnPTR).Adc[0] );
    
    //initialize y
    returnPTR -> y = malloc( points*sizeof(long double) );
    (*returnPTR).y[0] = powl( M*M_PI*freqInit/(300000000) , 1.0/3.0 );
    
    //calculate the initial values of the dynamic variables chieff, deltachi, and J
    long double LhatInit[3] , s1Init[3] , s2Init[3];
    LhatInit[0] = 0.0;
    LhatInit[1] = 0.0;
    LhatInit[2] = 1.0;
    s1Init[0] = (*returnPTR).S1*sinl( theta1 )/(*returnPTR).m1;
    s1Init[1] = 0.0;
    s1Init[2] = (*returnPTR).S1*cosl( theta1 )/(*returnPTR).m1;
    s2Init[2] = cosl( theta2 );
    s2Init[0] = ((*returnPTR).S2/(*returnPTR).m2)*( sinl(theta2)*cosl(theta12) );
    s2Init[1] = ((*returnPTR).S2/(*returnPTR).m2)*( sinl(theta2)*sinl(theta12) );
    s2Init[2] = ((*returnPTR).S2/(*returnPTR).m2)*( cosl(theta2) );
    long double deltachiInit = CalcDeltachi( s1Init , s2Init , LhatInit );
    long double chieffInit = CalcChieff( s1Init , s2Init , LhatInit );
    long double JInit = CalcJ( s1Init , s2Init , LhatInit , (*returnPTR).m1 , (*returnPTR).m2 , (*returnPTR).y[0] , (*returnPTR).eta );
    long double deltachiAMPSign = CalcdeltachiAMPSign( M , kappa1 , kappa2 , (*returnPTR).m1 , (*returnPTR).m2 , (*returnPTR).y[0] , s1Init , s2Init , LhatInit );
    long double chieffAMPSign = CalcchieffAMPSign( M , kappa1 , kappa2 , (*returnPTR).m1 , (*returnPTR).m2 , (*returnPTR).y[0] , s1Init , s2Init , LhatInit );
    returnPTR -> sign = deltachiAMPSign*chieffAMPSign;
    
    
    //initialize t
    long double yFinal , tFinal;
    yFinal = powl( M*M_PI*freqFinal/(300000000) , 1.0/3.0 );
    tFinal = (M*5.0/((*returnPTR).eta*256.0))*( powl( (*returnPTR).y[0] , -8.0 ) -  powl( yFinal , -8.0 ) );
    returnPTR -> t = malloc( points*sizeof(long double) );
    for( int i = 0; i < points; i++ )
    {
        (*returnPTR).t[i] = i*tFinal/(points-1);
    }
    printf("%Le  ,  %Le\n",(*returnPTR).y[0],yFinal);
    
    //initialize dcA, dcS, ceA, ceS, and psi
    returnPTR -> dcA = malloc( points*sizeof(long double) );
    returnPTR -> dcS = malloc( points*sizeof(long double) );
    returnPTR -> ceA = malloc( points*sizeof(long double) );
    returnPTR -> ceS = malloc( points*sizeof(long double) );
    returnPTR -> psi = malloc( points*sizeof(long double) );
    //precalculate the cubic coefficients to use in initial average/amplitude calculations (so we dont have to recaclulate in each function)
    long double B[2] , C[3] , D[4] , cubeRoots[3] , tempInit[4] , dydt , dpsidt , dJAVGdt;
    CalcB( JInit , (*returnPTR).eta/(*returnPTR).y[0] , (*returnPTR).S1 , (*returnPTR).S2 , (*returnPTR).y[0] , (*returnPTR).eta , (*returnPTR).dm , B );
    CalcC( JInit , (*returnPTR).eta/(*returnPTR).y[0] , (*returnPTR).S1 , (*returnPTR).S2 , (*returnPTR).y[0] , (*returnPTR).eta , (*returnPTR).dm , C );
    CalcD( JInit , (*returnPTR).eta/(*returnPTR).y[0] , (*returnPTR).S1 , (*returnPTR).S2 , (*returnPTR).y[0] , (*returnPTR).eta , (*returnPTR).dm , D );
    //calculate the initial averages and initial phase based on the first three derivatives of deltachi
    CalcInitAVGandpsi( (*returnPTR).Adc , (*returnPTR).Ace , B , C , D , deltachiInit , chieffInit , (*returnPTR).y[0] , (*returnPTR).dm , deltachiAMPSign , chieffAMPSign , (*returnPTR).dcA , (*returnPTR).ceA , (*returnPTR).psi );
    //calculate the amplitudes from the averages as well as the phase derivative and store them in the vector tempInit for retrieval
    CalcAMPsAndPhaseDeriv( (*returnPTR).Adc , (*returnPTR).Ace , B , C , D , (*returnPTR).dm , (*returnPTR).y[0] , (*returnPTR).dcA[0] , (*returnPTR).ceA[0] , deltachiAMPSign , chieffAMPSign , tempInit );
    (*returnPTR).dcS[0] = tempInit[0];
    (*returnPTR).ceS[0] = tempInit[1];
    dpsidt = tempInit[2];
    
    //calculate dydt for later
    dydt = Calcdydt( (*returnPTR).y[0] , (*returnPTR).eta );
    
    
    returnPTR -> JA = malloc( points*sizeof(long double) );
    returnPTR -> JS = malloc( points*sizeof(long double) );
    returnPTR -> JC = malloc( points*sizeof(long double) );
    //calculate the derivative of J average and the amplitudes
    CalcdJAVGAndJAMPs( dydt , dpsidt , (*returnPTR).eta/(*returnPTR).y[0] , (*returnPTR).y[0] , (*returnPTR).dm , JInit , (*returnPTR).dcA[0] , (*returnPTR).dcS[0] , (*returnPTR).ceA[0] , (*returnPTR).ceS[0] , tempInit );
    (*returnPTR).JS[0] = tempInit[1];
    (*returnPTR).JC[0] = tempInit[2];
    (*returnPTR).JA[0] = JInit - (*returnPTR).JS[0]*sinl((*returnPTR).psi[0]) - (*returnPTR).JC[0]*cosl((*returnPTR).psi[0]);
    dJAVGdt = tempInit[0];
    
    returnPTR -> phizA = malloc( points*sizeof(long double) );
    returnPTR -> phizC = malloc( points*sizeof(long double) );
    returnPTR -> zetaA = malloc( points*sizeof(long double) );
    returnPTR -> zetaC = malloc( points*sizeof(long double) );
    //evolve phiz and zeta
    CalcdphizAVGAndphizAMP( (*returnPTR).q , (*returnPTR).eta  , (*returnPTR).dm , (*returnPTR).S1 , (*returnPTR).S2 , (*returnPTR).y[0] , (*returnPTR).JA[0] , (*returnPTR).JS[0] , (*returnPTR).JC[0] , (*returnPTR).dcA[0] , (*returnPTR).dcS[0] , (*returnPTR).ceA[0] , (*returnPTR).ceS[0] , (*returnPTR).Ace , (*returnPTR).k1 , (*returnPTR).k2 , (*returnPTR).m1 , (*returnPTR).m2 , dpsidt , (*returnPTR).psi[0] , tempInit );
    (*returnPTR).phizA[0] = -tempInit[1];
    (*returnPTR).zetaA[0] = -tempInit[3];
    
    returnPTR -> phiz = malloc( points*sizeof(long double) );
    returnPTR -> zeta = malloc( points*sizeof(long double) );
    returnPTR -> thetaL = malloc( points*sizeof(long double) );
    (*returnPTR).phiz[0] = (*returnPTR).phizA[0] + tempInit[1];
    (*returnPTR).zeta[0] = (*returnPTR).zetaA[0] + tempInit[3];
    (*returnPTR).thetaL[0] = ANALCalcThetaL( JInit , (*returnPTR).eta/(*returnPTR).y[0] , (*returnPTR).ceA[0] + (*returnPTR).ceS[0]*sinl((*returnPTR).psi[0]) , (*returnPTR).dcA[0] + (*returnPTR).dcS[0]*sinl((*returnPTR).psi[0]) , (*returnPTR).dm );
    
    return returnPTR;
    
}


void freesys(  struct Asys* sys )
{
    free( (*sys).t );
    free( (*sys).y );
    free( (*sys).dcA );
    free( (*sys).dcS );
    free( (*sys).ceA );
    free( (*sys).ceS );
    free( (*sys).psi );
    free( (*sys).JA );
    free( (*sys).JS );
    free( (*sys).JC );
    free( (*sys).phizA );
    free( (*sys).phizC );
    free( (*sys).zetaA );
    free( (*sys).zetaC );
    free( (*sys).phiz );
    free( (*sys).zeta );
    free( (*sys).thetaL );
    free( sys );
}



//the return vector stores (in order) dydt, dpsidt, dJAVGdt, ddeltachiAVGdt, dchieffAVGdt, dphizAVGdt, dzetaAVGdt
void ANRK4STEP( struct Asys* sys , long double M , long double y , long double dcA , long double ceA , long double JA , long double returnVec[7] )
{
    //calculate dydt
    returnVec[0] = Calcdydt( y , (*sys).eta )/M;
    
    long double B[2] , C[3] , D[4] , cubeRoots[3] , tempInit[4] , dcS , ceS , JS , JC , phizC;
    CalcB( JA , (*sys).eta/y , (*sys).S1 , (*sys).S2 , y , (*sys).eta , (*sys).dm , B );
    CalcC( JA , (*sys).eta/y , (*sys).S1 , (*sys).S2 , y , (*sys).eta , (*sys).dm , C );
    CalcD( JA , (*sys).eta/y , (*sys).S1 , (*sys).S2 , y , (*sys).eta , (*sys).dm , D );
        
    //calculate dpsidt
    CalcAMPsAndPhaseDeriv( (*sys).Adc , (*sys).Ace , B , C , D , (*sys).dm , y , dcA , ceA , 1.0 , (*sys).sign , tempInit );
    dcS = tempInit[0];
    ceS = tempInit[1];
    returnVec[1] = tempInit[2]/M;
    
    //calculate dJAVGdt
    CalcdJAVGAndJAMPs( returnVec[0] , returnVec[1] , (*sys).eta/y , y , (*sys).dm , JA , dcA , dcS , ceA , ceS , tempInit  );
    JS = tempInit[1];
    JC = tempInit[2];
    returnVec[2] = tempInit[0];
    
    //calculate ddeltachiAVGdt, dchieffAVGdt
    CalcDeltachiAVGDeriv( B , C , D , (*sys).Adc , (*sys).Ace , (*sys).dm , (*sys).eta , (*sys).S1 , (*sys).S2 , y , returnVec[0] , JA , returnVec[2] , (*sys).eta/y , -(*sys).eta*returnVec[0]/(y*y) , dcA+dcS , dcA-dcS , ceA+ceS , ceA-ceS , tempInit );
    returnVec[3] = tempInit[0];
    returnVec[4] = tempInit[1];
    
    //calculate dphizAVGdt, dzetaAVGdt
    CalcdphizAVGAndphizAMP( (*sys).q , (*sys).eta  , (*sys).dm , (*sys).S1 , (*sys).S2 , y , JA , JS , JC , dcA , dcS , ceA , ceS , (*sys).Ace , (*sys).k1 , (*sys).k2 , (*sys).m1 , (*sys).m2 , 1.0 , 0.0 , tempInit );
    returnVec[5] = tempInit[0]/M;
    returnVec[6] = tempInit[2]/M;
    
}



void RK4derivsum( long double k1[7] , long double k2[7] , long double k3[7] , long double k4[7] , long double k[7] )
{
    for( int i = 0; i < 7; i++ )
    {
        k[i] = ( k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i] )/6.0;
    }
}



//Function to semi-analytically evolve the system
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
struct Waveform* AnalyticalEvolver( long double freqInit , long double freqFinal , long double M , long double q , long double chi1 , long double chi2 , long double theta1 , long double theta2 , long double theta12 , long double kappa1 , long double kappa2 , int points )
{
    
    struct Asys* sys = Initsys( points , M , freqInit , freqFinal , q , chi1 , chi2 , theta1 , theta2 , theta12 , kappa1 , kappa2 );

    long double k1[7] , k2[7] , k3[7] , k4[7] , k[7];
    long double B[2] , C[3] , D[4] , cubeRoots[3] , tempInit[4] , dydt , dpsidt;
    for( int i = 1; i < points; i++ )
    {
        ANRK4STEP( sys , M , (*sys).y[i-1] , (*sys).dcA[i-1] , (*sys).ceA[i-1] , (*sys).JA[i-1] , k1 );
        ANRK4STEP( sys , M , (*sys).y[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*0.5*k1[0] , (*sys).dcA[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*0.5*k1[3] , (*sys).ceA[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*0.5*k1[4] , (*sys).JA[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*0.5*k1[2] , k2 );
        ANRK4STEP( sys , M , (*sys).y[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*0.5*k2[0] , (*sys).dcA[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*0.5*k2[3] , (*sys).ceA[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*0.5*k2[4] , (*sys).JA[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*0.5*k2[2] , k3 );
        ANRK4STEP( sys , M , (*sys).y[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*k3[0] , (*sys).dcA[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*k3[3] , (*sys).ceA[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*k3[4] , (*sys).JA[i-1] + ( (*sys).t[i] - (*sys).t[i-1] )*k3[2] , k4 );
        RK4derivsum( k1 , k2 , k3 , k4 , k );
        
        //update the averages using results of Rk4
        (*sys).y[i] = k[0]*( (*sys).t[i] - (*sys).t[i-1] ) + (*sys).y[i-1];
        (*sys).psi[i] = k[1]*( (*sys).t[i] - (*sys).t[i-1] ) + (*sys).psi[i-1];
        (*sys).JA[i] = k[2]*( (*sys).t[i] - (*sys).t[i-1] ) + (*sys).JA[i-1];
        (*sys).dcA[i] = k[3]*( (*sys).t[i] - (*sys).t[i-1] ) + (*sys).dcA[i-1];
        (*sys).ceA[i] = k[4]*( (*sys).t[i] - (*sys).t[i-1] ) + (*sys).ceA[i-1];
        (*sys).phizA[i] = k[5]*( (*sys).t[i] - (*sys).t[i-1] ) + (*sys).phizA[i-1];
        (*sys).zetaA[i] = k[6]*( (*sys).t[i] - (*sys).t[i-1] ) + (*sys).zetaA[i-1];
        
        //find the new amplitudes from the new averages
        //start by calculating current coefficients
        CalcB( (*sys).JA[i] , (*sys).eta/(*sys).y[i] , (*sys).S1 , (*sys).S2 , (*sys).y[i] , (*sys).eta , (*sys).dm , B );
        CalcC( (*sys).JA[i] , (*sys).eta/(*sys).y[i] , (*sys).S1 , (*sys).S2 , (*sys).y[i] , (*sys).eta , (*sys).dm , C );
        CalcD( (*sys).JA[i] , (*sys).eta/(*sys).y[i] , (*sys).S1 , (*sys).S2 , (*sys).y[i] , (*sys).eta , (*sys).dm , D );
        
        //then calculate amps and current phase deriv for dc and ce
        CalcAMPsAndPhaseDeriv( (*sys).Adc , (*sys).Ace , B , C , D , (*sys).dm , (*sys).y[i] , (*sys).dcA[i] , (*sys).ceA[i] , 1.0 , (*sys).sign , tempInit );
        (*sys).dcS[i] = tempInit[0];
        (*sys).ceS[i] = tempInit[1];
        dpsidt = tempInit[2];
        dydt = Calcdydt( (*sys).y[i] , (*sys).eta );
        
        //calculate the amps for J
        CalcdJAVGAndJAMPs( dydt , dpsidt , (*sys).eta/(*sys).y[i] , (*sys).y[i] , (*sys).dm , (*sys).JA[i] , (*sys).dcA[i] , (*sys).dcS[i] , (*sys).ceA[i] , (*sys).ceS[i] , tempInit  );
        (*sys).JS[i] = tempInit[1];
        (*sys).JC[i] = tempInit[2];
        
        //calculate the amps for phiz and zeta
        CalcdphizAVGAndphizAMP( (*sys).q , (*sys).eta  , (*sys).dm , (*sys).S1 , (*sys).S2 , (*sys).y[i] , (*sys).JA[i] , (*sys).JS[i] , (*sys).JC[i] , (*sys).dcA[i] , (*sys).dcS[i] , (*sys).ceA[i] , (*sys).ceS[i] , (*sys).Ace , (*sys).k1 , (*sys).k2 , (*sys).m1 , (*sys).m2 , dpsidt , (*sys).psi[i] , tempInit );
        
        //calculate zeta, phiz, and thetaL for the current point.
        (*sys).phiz[i] = (*sys).phizA[i] + tempInit[1];
        (*sys).zeta[i] = (*sys).zetaA[i] + tempInit[3];
        (*sys).thetaL[i] = ANALCalcThetaL( (*sys).JA[i] + (*sys).JS[i]*sinl( (*sys).psi[i] ) + (*sys).JC[i]*cosl( (*sys).psi[i] ) , (*sys).eta/(*sys).y[i] , (*sys).ceA[i] + (*sys).ceS[i]*sinl((*sys).psi[i]) , (*sys).dcA[i] + (*sys).dcS[i]*sinl((*sys).psi[i]) , (*sys).dm );

    }
    
    
    //output the relevant data in csv form to a text document for plotting purposes
    char filename[100];
    sprintf(filename,"plotAnalytical - k1 = %.2Lf , k2 = %.2Lf.txt" , kappa1 , kappa2);
    FILE *filey = fopen(filename, "w");
    fprintf(filey, "value number,time,chieff,deltachi,J,phiz,thetaL \n");
    for( int i = 0; i < points; i++ )
    {
        fprintf(filey, "val#%i,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le,%.9Le \n", i ,  (*sys).t[i], (*sys).ceA[i]+(*sys).ceS[i]*sinl((*sys).psi[i]) , (*sys).dcA[i]+(*sys).dcS[i]*sinl((*sys).psi[i])  , (*sys).JA[i] + (*sys).JS[i]*sinl((*sys).psi[i]) + (*sys).JC[i]*cosl((*sys).psi[i]) , (*sys).phiz[i] , (*sys).zeta[i] , (*sys).thetaL[i] , (*sys).dcA[i] , (*sys).ceA[i] );

    }
    fflush(filey);
    
    
    struct Waveform* returner = GenWaveform( M , (*sys).t , (*sys).y , (*sys).phiz , (*sys).zeta , (*sys).thetaL , points , 0 , 0 , 0 , 0 , (*sys).eta );
    freesys( sys );
    return returner;

}
