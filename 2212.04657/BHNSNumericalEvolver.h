//
//  BHNSNumericalEvolver.h
//  
//
//  Created by Michael Lahaye on 2022-04-02.
//

#ifndef BHNSNumericalEvolver_h
#define BHNSNumericalEvolver_h

#include <stdio.h>

extern long double Calcdydt( long double y , long double eta );

extern long double CalcDeltachi( long double s1[3] , long double s2[3] , long double Lhat[3] );
extern long double CalcChieff( long double s1[3] , long double s2[3] , long double Lhat[3] );
extern long double CalcJ( long double s1[3] , long double s2[3] , long double Lhat[3] , long double m1 , long double m2 , long double y , long double eta );

extern void CalcAngularMomentaDerivs( long double M , long double kappa1 , long double kappa2 , long double m1 , long double m2 , long double y , long double s1[3] , long double s2[3] , long double Lhat[3] , long double returnVectors1[3] , long double returnVectors2[3] , long double returnVectorLhat[3] );

extern struct Waveform* NumericalEvolver( long double freqInit , long double freqFinal , long double M , long double q , long double chi1 , long double chi2 , long double theta1 , long double theta2 , long double theta12 , long double kappa1 , long double kappa2 , int points );

#endif /* BHNSNumericalEvolver_h */
