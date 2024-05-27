//------------------------------------------------------------------------------
//
// Geodetic (class implementation)
//
// Purpose:
//
//   Class (with all public elements) for handling geodetic coordinates
//
//------------------------------------------------------------------------------

// Default constructor

Geodetic::Geodetic ()
 : lon(0.0), lat(0.0), h(0.0)
{
}

// Simple constructor

Geodetic::Geodetic (double lambda, double phi, double alt)
{
  lon=lambda; lat=phi; h=alt;
}

// Constructor for geodetic coordinates from given position

Geodetic::Geodetic (Vector r,                         // Position vector [m]
                    double R_equ,                     // Equator radius [m]
                    double f)                         // Flattening
{

  const double  eps     = 1.0e3*eps_mach;   // Convergence criterion 
  const double  epsRequ = eps*R_equ;
  const double  e2      = f*(2.0-f);        // Square of eccentricity
  
  const double  X = r(0);                   // Cartesian coordinates
  const double  Y = r(1);
  const double  Z = r(2);
  const double  rho2 = X*X + Y*Y;           // Square of distance from z-axis
  
  // Check validity of input data
  
  if (Norm(r)==0.0) {
    cerr << " invalid input in Geodetic constructor" << endl;
    lon=0.0; lat=0.0; h=-R_Earth;
    return;
  }

  // Iteration 

  double  dZ, dZ_new, SinPhi;
  double  ZdZ, Nh, N;

  dZ = e2*Z;
  for(;;) {
    ZdZ    =  Z + dZ;
    Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
    SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
    N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi); 
    dZ_new =  N*e2*SinPhi;
    if ( fabs(dZ-dZ_new) < epsRequ ) break;
    dZ = dZ_new;
  }
    
  // Longitude, latitude, altitude

  lon = atan2 ( Y, X );
  lat = atan2 ( ZdZ, sqrt(rho2) );
  h   = Nh - N;

}

// Position vector [m] from geodetic coordinates

Vector Geodetic::Position (double R_equ,   // Equator radius [m]
                           double f    )   // Flattening
  const
{  

  const double  e2     = f*(2.0-f);        // Square of eccentricity
  const double  CosLat = cos(lat);         // (Co)sine of geodetic latitude
  const double  SinLat = sin(lat);

  double  N;
  Vector  r(3);
      
  // Position vector 

  N = R_equ / sqrt(1.0-e2*SinLat*SinLat);

  r(0) =  (         N+h)*CosLat*cos(lon);
  r(1) =  (         N+h)*CosLat*sin(lon);
  r(2) =  ((1.0-e2)*N+h)*SinLat;

  return r;

}

// Transformation to local tangent coordinates

Matrix Geodetic::LTC_Matrix () const
{
  return LTCMatrix(lon,lat);
}
