
//------------------------------------------------------------------------------
//
// AzEl
//
// Purpose:
//
//   Computes azimuth and elevation from local tangent coordinates
//
// Input/Output:
//
//   s   Topocentric local tangent coordinates (East-North-Zenith frame)
//   A   Azimuth [rad]
//   E   Elevation [rad]
//
//------------------------------------------------------------------------------

void AzEl ( const Vector& s, double& A, double& E ) 
{
  A = atan2(s(0),s(1));
  A = ((A<0.0)? A+pi2 : A);
  E = atan ( s(2) / sqrt(s(0)*s(0)+s(1)*s(1)) );
}


//------------------------------------------------------------------------------
//
// AzEl
//
// Purpose:
//
//   Computes azimuth, elevation and partials from local tangent coordinates
//
// Input/Output:
//
//   s      Topocentric local tangent coordinates (East-North-Zenith frame)
//   A      Azimuth [rad]
//   E      Elevation [rad]
//   dAds   Partials of azimuth w.r.t. s
//   dEds   Partials of elevation w.r.t. s
//
//------------------------------------------------------------------------------

void AzEl( const Vector& s, double& A, double& E, Vector& dAds, Vector& dEds ) 
{
  const double rho = sqrt(s(0)*s(0)+s(1)*s(1));
  // Angles
  A = atan2(s(0),s(1));
  A = ((A<0.0)? A+pi2 : A);
  E = atan ( s(2) / rho );
  // Partials
  dAds = Vector( s(1)/(rho*rho), -s(0)/(rho*rho), 0.0 );
  dEds = Vector( -s(0)*s(2)/rho, -s(1)*s(2)/rho , rho ) / Dot(s,s);
}
