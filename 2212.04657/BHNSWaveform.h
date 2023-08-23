//
//  BHNSWaveform.h
//  
//
//  Created by Michael Lahaye on 2022-04-17.
//

#ifndef BHNSWaveform_h
#define BHNSWaveform_h

#include <complex.h>
#include <stdio.h>


struct Waveform
{
    long double *f;
    long double complex *h;
};

struct Waveform* GenWaveform( long double M , long double *tValues , long double *y , long double *phiz , long double *thetaL , long double *zeta , int points , long double thetaS , long double phiS , long double thetaN , long double phiN , long double eta );

void freeWaveform( struct Waveform* PTR );

long double CalcMismatch( struct Waveform* Wvfrm1 , struct Waveform* Wvfrm2 , int points , int ratio );

#endif /* BHNSWaveform_h */
