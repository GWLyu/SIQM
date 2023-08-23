//
//  BHNSMain.c
//
//
//  Created by Michael Lahaye on 2022-04-02.
//

#include "BHNSMain.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "BHNSNumericalEvolver.h"
#include "BHNSAnalyticalEvolver.h"
#include "BHNSWaveform.h"
#include <complex.h>

//the main
int main(){

    int points = 10000;
    int ratio = 20;
    long double fInit = 10;
    long double fFinal = 1.0*powl(10.0,2.0);
    long double M = 1476.625;
    long double q = 0.5;
    struct Waveform *NUMNSNS , *NUMBHBH , *ANNSNS , *ANBHBH;
    long double NUMBHBHvsNUMNSNS, ANNSNSvsNUMNSNS;
    int i = 1;
    
    FILE *filey = fopen("mismatch.txt", "w");
    fprintf(filey, " ,q,χ1,χ2,κ1,κ2,θ1,θ2,θ3,mismatch vs NUM BHBH, mismatch vs AN NSNS \n");
    
        for( long double chi1 = 0.3; chi1 < 0.4; chi1+=0.1 )
        {
            for( long double chi2 = 0.3; chi2 < 0.4; chi2+=0.1 )
            {
                for( long double kappa1 = 1.0; kappa1 < 1.2; kappa1+=0.2 )
                {
                    for( long double kappa2 = 6.0; kappa2 < 6.5; kappa2+=0.2 )
                    {
                        for( long double theta1coeff = 0.1; theta1coeff < 1.01; theta1coeff+=0.5 )
                        {
                            for( long double theta2coeff = 0.1; theta2coeff < 1.01; theta2coeff+=0.5 )
                            {
                                for( long double theta3coeff = 0.1; theta3coeff < 2.0; theta3coeff+=0.8 )
                                {
                                    NUMNSNS = NumericalEvolver( fInit , fFinal , M, q , chi1 , chi2 , theta1coeff*M_PI , theta2coeff*M_PI , theta3coeff*M_PI , kappa1 , kappa2 , (points-1)*ratio+1 );
                                    NUMBHBH = NumericalEvolver( fInit , fFinal , M, q , chi1 , chi2 , theta1coeff*M_PI , theta2coeff*M_PI , theta3coeff*M_PI , 1.0 , 1.0 , (points-1)*ratio+1 );
                                    ANNSNS = AnalyticalEvolver( fInit , fFinal , M , q , chi1 , chi2 , theta1coeff*M_PI , theta2coeff*M_PI , theta3coeff*M_PI , kappa1 , kappa2 , points );
                                    //ANBHBH = AnalyticalEvolver( fInit , fFinal , 1.0 , q , chi1 , chi2 , theta1coeff*M_PI , theta2coeff*M_PI , theta3coeff*M_PI , 1.0 , 1.0 , points );
                                    
                                    NUMBHBHvsNUMNSNS = CalcMismatch( NUMNSNS , NUMBHBH , (points-1)*ratio+1 , 1 );
                                    if( NUMBHBHvsNUMNSNS != NUMBHBHvsNUMNSNS )
                                        NUMBHBHvsNUMNSNS = 0.0;
                                    ANNSNSvsNUMNSNS = CalcMismatch( ANNSNS , NUMNSNS , points , ratio );
                                    if( ANNSNSvsNUMNSNS != ANNSNSvsNUMNSNS )
                                        ANNSNSvsNUMNSNS = 0.0;
                                    fprintf(filey, "case#%i,%.8Le,%.8Le,%.8Le,%.8Le,%.8Le,%.8Le,%.8Le,%.8Le,%.8Le,%.8Le \n", i , q , chi1 , chi2 , kappa1 , kappa2 , theta1coeff*M_PI , theta2coeff*M_PI , theta3coeff*M_PI , NUMBHBHvsNUMNSNS , ANNSNSvsNUMNSNS );
                                    
                                    freeWaveform( NUMNSNS );
                                    freeWaveform( NUMBHBH );
                                    freeWaveform( ANNSNS );
                                    //freeWaveform( ANBHBH );
                                    
                                    printf( "\rDone waveform %d    ", i );
                                    fflush(stdout);
                                    i++;
                                }
                            }
                        }
                    }
                }
            }
        }
  
    fflush(filey);
    printf( "\rProcess Complete   \n" );
    
    /*
    FILE *filey = fopen("frequencies.txt", "w");
    fprintf(filey, " ,f,|hNSNS/hBHBH|,arg(hNSNS/hBHBH) \n");
    long double deltachi , chieff , J;
    for( int i = 0; i < points; i++ )
    {
        fprintf(filey, "val#%i,%.8Le,%.8Le,%.8Le \n", i , ((*ANNSNS).f)[i] , cabsl( ((*ANNSNS).h)[i]/((*NUMNSNS).h)[i] ) , cargl( ((*ANNSNS).h)[i]/((*NUMNSNS).h)[i] ) );
    }
    fflush(filey);
    */
    
    return 0;
}
