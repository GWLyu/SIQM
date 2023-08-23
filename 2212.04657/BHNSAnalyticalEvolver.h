//
//  BHNSAnalyticEvolver.h
//  
//
//  Created by Michael Lahaye on 2022-04-02.
//

#ifndef BHNSAnalyticalEvolver_h
#define BHNSAnalyticalEvolver_h

#include <stdio.h>

extern struct Waveform* AnalyticalEvolver( long double freqInit , long double freqFinal , long double M , long double q , long double chi1 , long double chi2 , long double theta1 , long double theta2 , long double theta12 , long double kappa1 , long double kappa2 , int points );

struct Asys
{
    long double M , q , m1 , m2 , dm , eta , S1 , S2 , k1 , k2 , sign;
    long double Adc[2] , Ace[2];
    long double *t , *y , *psi , *dcA , *dcS , *ceA , *ceS , *JA , *JS , *JC , *phizA , *phizC , *phiz , *thetaL , *zetaA , *zetaC , *zeta ;
};

#endif /* BHNSAnalyticEvolver_h */
