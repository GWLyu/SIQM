//
//  BHNSVector.h
//
//
//  Created by Michael Lahaye on 2022-02-12.
//

#ifndef BHNSVectors_h
#define BHNSVectors_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



extern void ScalarMultiply( long double scalar , long double vector[3] , long double returnVector[3] );

extern void SumVectors( long double vector1[3] , long double vector2[3] , long double returnVector[3] );

extern long double DotProduct( long double vector1[3] , long double vector2[3] );

extern void CrossProduct( long double vector1[3] , long double vector2[3] , long double returnVector[3] );

extern long double Magnitude( long double vector[3] );

extern void Normalize( long double vector[3] );

extern void MatrixMultiply( long double matrix[3][3] , long double vector[3] , long double returnVector[3] );



#endif /* BHNSVector_h */
