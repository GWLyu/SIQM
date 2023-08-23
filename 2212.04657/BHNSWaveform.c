//
//  BHNSWaveform.c
//
//
//  Created by Michael Lahaye on 2022-04-17.
//

#include "BHNSWaveform.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>


long double complex WignerDP( int mprime , long double phiz , long double thetaL , long double zeta )
{
    long double complex temp;
    switch( mprime )
    {
        case -2:
            temp = 0.25*(1.0-cosl(thetaL))*(1.0-cosl(thetaL));
            break;
        case -1:
            temp = -0.5*sinl(thetaL)*(1.0+cosl(thetaL));
            break;
        case 0:
            temp = powl(3/8,0.5)*sinl(thetaL)*sinl(thetaL);
            break;
        case 1:
            temp = -0.5*sinl(thetaL)*(1.0-cosl(thetaL));
            break;
        case 2:
            temp = 0.25*(1.0+cosl(thetaL))*(1.0+cosl(thetaL));
            break;
    }
    return cexpl( -I*mprime*phiz )*temp*cexpl( -I*2.0*zeta );
}



long double complex WignerDM( int mprime , long double phiz , long double thetaL , long double zeta )
{
    long double complex temp;
    switch( mprime )
    {
        case -2:
            temp = 0.25*(1.0+cosl(thetaL))*(1.0+cosl(thetaL));
            break;
        case -1:
            temp = -0.5*sinl(thetaL)*(1.0+cosl(thetaL));
            break;
        case 0:
            temp = powl(3/8,0.5)*sinl(thetaL)*sinl(thetaL);
            break;
        case 1:
            temp = -0.5*sinl(thetaL)*(1.0-cosl(thetaL));
            break;
        case 2:
            temp = 0.25*(1.0-cosl(thetaL))*(1.0-cosl(thetaL));
            break;
    }
    return cexpl( -I*mprime*phiz )*temp*cexpl( I*2.0*zeta );
}



long double complex SWSH( int mprime , long double thetaS , long double phiS )
{
    long double complex temp;
    switch( mprime )
    {
        case -2:
            temp = 0.25*(1.0-cosl(thetaS))*(1.0-cosl(thetaS));
            break;
        case -1:
            temp = -0.5*sinl(thetaS)*(1.0-cosl(thetaS));
            break;
        case 0:
            temp = powl(3/8,0.5)*sinl(thetaS)*sinl(thetaS);
            break;
        case 1:
            temp = -0.5*sinl(thetaS)*(1.0+cosl(thetaS));
            break;
        case 2:
            temp = 0.25*(1.0+cosl(thetaS))*(1.0+cosl(thetaS));
            break;
    };
    return powl(-1.0,mprime)*powl( 4*M_PI/5 , -0.5 )*cexpl( I*mprime*phiS )*temp;
}



long double complex hlm( int m , long double eta , long double y )
{
    long double complex temp;
    long double e;
    e = 0.57721566490153286060;
    if(m==2)
    {
        temp = 1;
        temp += y*( -107.0/42.0 + eta*55.0/42.0 );
        temp += 2*M_PI*powl(y,1.5);
        temp += y*y*(-2173.0/1512.0 - eta*1069.0/216.0 + eta*eta*2047.0/1512.0);
        temp += powl(y,2.5)*( -107.0*M_PI/21.0 - 24.0*I*eta + eta*34.0*M_PI/21.0 );
    }
    else
    {
        temp = 1;
        temp += y*( -107.0/42.0 + eta*55.0/42.0 );
        temp += 2*M_PI*powl(y,1.5);
        temp += y*y*(-2173.0/1512.0 - eta*1069.0/216.0 + eta*eta*2047.0/1512.0);
        temp += powl(y,2.5)*( -107.0*M_PI/21.0 + 24.0*I*eta + eta*34.0*M_PI/21.0 );
    }
    return temp;
}



struct Waveform* GenWaveform( long double M , long double *tValues , long double *y , long double *phiz , long double *thetaL , long double *zeta , int points , long double thetaS , long double phiS , long double thetaN , long double phiN , long double eta )
{

    struct Waveform* returnPTR = malloc( sizeof(struct Waveform) );
    returnPTR -> f = malloc( points*sizeof(long double) );
    returnPTR -> h = malloc( points*sizeof(long double complex) );
    
    
    long double Phi , phiOrb , Tm , psip;
    long double Fplus , Fcross;
    long double complex Hlm , temp1 , temp2;
    phiOrb = 0;
    Phi = phiOrb - 3.0*y[0]*y[0]*y[0]*(2.0-eta*y[0]*y[0])*logl( y[0] );
    
    psip = 0;
    Fplus = 0.5*( 1.0 + cosl(thetaN)*cosl(thetaN) )*cosl(2.0*phiN)*cosl(2*psip) - cosl(thetaN)*sinl(2.0*phiN)*sinl(2*psip);
    Fcross = 0.5*( 1.0 + cosl(thetaN)*cosl(thetaN) )*cosl(2.0*phiN)*sinl(2*psip) - cosl(thetaN)*sinl(2.0*phiN)*cosl(2*psip);
    
    
    for( int i = 0; i < points; i++ )
    {
        if(i!=0)
        {
            phiOrb = phiOrb + (tValues[i]-tValues[i-1])*y[i-1]*y[i-1]*y[i-1];
            Phi = phiOrb - 3.0*y[i]*y[i]*y[i]*(2.0-eta*y[i]*y[i])*logl( y[i] );
            
            Tm = powl( 6.0*y[i]*y[i]*(y[i]-y[i-1])/(tValues[i]-tValues[i-1]) , -0.5 );
        }
        else
        {
            Tm = powl( 6.0*y[i+1]*y[i+1]*(y[i+1]-y[i])/(tValues[i+1]-tValues[i]) , -0.5 );
        }
        
        
        ((*returnPTR).f)[i] = (300000000)*powl( y[i] , 3.0 )/(M*M_PI);
        
        temp1 = SWSH( -2 , thetaS , phiS )*WignerDP( -2 , phiz[i] , thetaL[i] , zeta[i] ) + SWSH( -1 , thetaS , phiS )*WignerDP( -1 , phiz[i] , thetaL[i] , zeta[i] ) + SWSH( 0 , thetaS , phiS )*WignerDP( 0 , phiz[i] , thetaL[i] , zeta[i] ) + SWSH( 1 , thetaS , phiS )*WignerDP( 1 , phiz[i] , thetaL[i] , zeta[i] ) + SWSH( 2 , thetaS , phiS )*WignerDP( 2 , phiz[i] , thetaL[i] , zeta[i] );
        temp2 = SWSH( -2 , thetaS , phiS )*WignerDM( -2 , phiz[i] , thetaL[i] , zeta[i] ) + SWSH( -1 , thetaS , phiS )*WignerDM( 1 , phiz[i] , thetaL[i] , zeta[i] ) + SWSH( 0 , thetaS , phiS )*WignerDM( 0 , phiz[i] , thetaL[i] , zeta[i] ) + SWSH( 1 , thetaS , phiS )*WignerDM( 1 , phiz[i] , thetaL[i] , zeta[i] ) + SWSH( 2 , thetaS , phiS )*WignerDM( 2 , phiz[i] , thetaL[i] , zeta[i] );
        
        Hlm = 0.5*( Fplus + I*Fcross )*temp1*hlm( 2 , eta , y[i] ) + 0.5*( Fplus - I*Fcross )*temp2*hlm( -2 , eta , y[i] );
        
        ((*returnPTR).h)[i] = Hlm*powl( 2.0*M_PI , 0.5 )*Tm;//cexpl( I*( 2.0*M_PI*((*returnPTR).f)[i]*tValues[i] - 2.0*Phi - M_PI/4.0 ) );
    }
    
    /*
    for( int i = 0; i < points; i++ )
    {
        ((*returnPTR).f)[i] = powl( y[i] , 3.0 )/M_PI;
        ((*returnPTR).h)[i] = cexpl( I*( phiz[i] ) );
    }
    */
    return returnPTR;
    
}


void freeWaveform( struct Waveform* PTR )
{
    free( (*PTR).h );
    free( (*PTR).f );
    free( PTR );
}


long double CalcInnerProd( struct Waveform* Wvfrm1 , struct Waveform* Wvfrm2 , int points , int ratio )
{
    int j;
    long double complex result = 0.0;
    for( int i = 1; i < points; i++)
    {
        j = ratio*i;
        result += 0.5*( ((*Wvfrm1).h)[i]*conjl( ((*Wvfrm2).h)[j] ) + ((*Wvfrm1).h)[i-1]*conjl( ((*Wvfrm2).h)[j-ratio] ) )*( ((*Wvfrm1).f)[i] - ((*Wvfrm1).f)[i-1] );
    }
    return 4.0*creall( result );
}


long double CalcMismatch( struct Waveform* Wvfrm1 , struct Waveform* Wvfrm2 , int points , int ratio )
{
    long double A , B , C;
    A = CalcInnerProd( Wvfrm1 , Wvfrm2 , points , ratio );
    B = CalcInnerProd( Wvfrm1 , Wvfrm1 , points , 1 );
    C = CalcInnerProd( Wvfrm2 , Wvfrm2 , (points-1)*ratio+1 , 1 );
    return A/powl( B*C , 0.5 );
}
