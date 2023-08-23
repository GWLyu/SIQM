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
    
    long double kappa1 = 1.0;
    long double kappa2 = 2.5;
    int points = 10000;
    int ratio = 20;
    long double fInit = 10;
    long double fFinal = 1.0*powl(10.0,3.0);
    //long double M = (2.0+2.6)*2900.0;
    long double m1 = 2.6;
    long double m2 = 1.5;
    long double M = (m1 + m2)*1476.625;
    long double q = m2/m1;
    long double theta1 = 0.05*M_PI;
    long double theta2 = 0.25*M_PI;
    long double theta12 = 0.1*M_PI;
    struct Waveform *NUMNSNS , *NUMBHBH , *ANNSNS , *ANBHBH;
    NUMNSNS = NumericalEvolver( fInit , fFinal , M , q , 0.4 , 0.7 , theta1 , theta2 , theta12 , kappa1 , kappa2 , ratio*(points-1)+1 );
    NUMBHBH = NumericalEvolver( fInit , fFinal , M , q , 0.4 , 0.7 , theta1 , theta2 , theta12 , 1.0 , 1.0 , ratio*(points-1)+1 );
    ANNSNS = AnalyticalEvolver( fInit , fFinal , M , q , 0.4 , 0.7 , theta1 , theta2 , theta12 , kappa1 , kappa2 , points );
    ANBHBH = AnalyticalEvolver( fInit , fFinal , M , q , 0.4 , 0.7 , theta1 , theta2 , theta12 , 1.0 , 1.0 , points );
    //0.5,0.4,0.5 the one that illustrates the very weird property
    //0.2,0.4,0.5
    //0.5,0.8,0.5 the first one
    //0.1,0.9,0.85 the physically realistic one
    
    FILE *filey = fopen("frequencies1.txt", "w");
    fprintf(filey, " ,f,|hNSNS/hBHBH|,arg(hNSNS/hBHBH) \n");
    for( int i = 0; i < points; i++ )
    {
        fprintf(filey, "val#%i,%.9Le,%.9Le,%.9Le \n", i , ((*ANNSNS).f)[i] , cabsl( ((*ANNSNS).h)[i]/((*NUMNSNS).h)[ratio*i] ) , cargl( ((*ANNSNS).h)[i]/((*NUMNSNS).h)[ratio*i] ) );
    }
    fflush(filey);
    
    FILE *filey2 = fopen("frequencies2.txt", "w");
    fprintf(filey2, " ,f,|hNSNS/hBHBH|,arg(hNSNS/hBHBH) \n");
    for( int i = 0; i < ratio*(points-1)+1; i++ )
    {
        fprintf(filey2, "val#%i,%.9Le \n", i , ((*NUMNSNS).f)[i] );
    }
    fflush(filey);
    
    return 0;
}
