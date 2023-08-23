//
//  BHNSVector.c
//
//
//  Created by Michael Lahaye on 2022-02-12.
//

#include "BHNSVectors.h"



//Calculates the scalar multiple of a vector.
//Inputs:
//      -the scalar
//      -the vector
//      -location to return the result
//Outputs:
//      -none
void ScalarMultiply( long double scalar , long double vector[3] , long double returnVector[3] )
{
    returnVector[0] = scalar*vector[0];
    returnVector[1] = scalar*vector[1];
    returnVector[2] = scalar*vector[2];
}



//Calculates the sum of two vectors.
//Inputs:
//      -the first vector
//      -the second vector
//      -location to return the result
//Outputs:
//      -none
void SumVectors( long double vector1[3] , long double vector2[3] , long double returnVector[3] )
{
    returnVector[0] = vector1[0] + vector2[0];
    returnVector[1] = vector1[1] + vector2[1];
    returnVector[2] = vector1[2] + vector2[2];
}



//Calculates the dot product of two vectors.
//Inputs:
//      -the first vector
//      -the second vector
//Outputs:
//      -the dot product
long double DotProduct( long double vector1[3] , long double vector2[3] )
{
    return vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2];
}



//Calculates the cross product of two vectors.
//Inputs:
//      -the first vector
//      -the second vector
//      -location to return the result
//Outputs:
//      -none
void CrossProduct( long double vector1[3] , long double vector2[3] , long double returnVector[3] )
{
    long double temp1 , temp2 , temp3;
    temp1 = vector1[1]*vector2[2] - vector1[2]*vector2[1];
    temp2 = vector1[2]*vector2[0] - vector1[0]*vector2[2];
    temp3 = vector1[0]*vector2[1] - vector1[1]*vector2[0];
    returnVector[0] = temp1;
    returnVector[1] = temp2;
    returnVector[2] = temp3;
}



//Calculates the magnitude of a vector.
//Inputs:
//      -the vector
//Outputs:
//      -the magnitude of the vector
long double Magnitude( long double vector[3] )
{
    return powl( ( vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2] ) , 0.5 );
}



//Normalizes a vector in place.
//Inputs:
//      -the vector
//Outputs:
//      -none
void Normalize( long double vector[3] )
{
    long double magnitude = Magnitude( vector );
    vector[0] = vector[0]/magnitude;
    vector[1] = vector[1]/magnitude;
    vector[2] = vector[2]/magnitude;
}



//Multiplies a 3x3 matrix by a vector.
//Inputs:
//      -the matrix
//      -the vector
//      -location to return the result
//Outputs:
//      -none
void MatrixMultiply( long double matrix[3][3] , long double vector[3] , long double returnVector[3] )
{
    returnVector[0] = matrix[0][0]*vector[0] + matrix[0][1]*vector[1] + matrix[0][2]*vector[2];
    returnVector[1] = matrix[1][0]*vector[0] + matrix[1][1]*vector[1] + matrix[1][2]*vector[2];
    returnVector[2] = matrix[2][0]*vector[0] + matrix[2][1]*vector[1] + matrix[2][2]*vector[2];
}
