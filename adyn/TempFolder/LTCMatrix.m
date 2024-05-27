// Transformation to local tangent coordinates

Matrix Geodetic::LTC_Matrix () const
{
  return LTCMatrix(lon,lat);
}

//------------------------------------------------------------------------------
//
// LTCMatrix
//
// Purpose:
//
//   Transformation from Greenwich meridian system to local tangent coordinates
//
// Input/Output:
//
//   lambda    Geodetic East longitude [rad]
//   phi       Geodetic latitude [rad]
//   <return>  Rotation matrix from the Earth equator and Greenwich meridian
//             to the local tangent (East-North-Zenith) coordinate system
//
//------------------------------------------------------------------------------

Matrix LTCMatrix (double lambda, double phi)
{
  
  Matrix  M(3,3);
  double  Aux;
  
  // Transformation to Zenith-East-North System
  M = R_y(-phi)*R_z(lambda);
  
  // Cyclic shift of rows 0,1,2 to 1,2,0 to obtain East-North-Zenith system
  for (int j=0; j<3; j++) {
    Aux=M(0,j); M(0,j)=M(1,j); M(1,j)=M(2,j); M(2,j)= Aux;
  }
  
  return  M;
  
}

