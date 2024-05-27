//------------------------------------------------------------------------------
//
// GMST
//
// Purpose:
//
//   Greenwich Mean Sidereal Time
//
// Input/Output:
//
//   Mjd_UT1   Modified Julian Date UT1
//   <return>  GMST in [rad]
//
//------------------------------------------------------------------------------

double GMST (double Mjd_UT1)
{

  // Constants

  const double Secs = 86400.0;        // Seconds per day

  // Variables

  double Mjd_0,UT1,T_0,T,gmst;

  // Mean Sidereal Time
  
  Mjd_0 = floor(Mjd_UT1);
  UT1   = Secs*(Mjd_UT1-Mjd_0);          // [s]
  T_0   = (Mjd_0  -MJD_J2000)/36525.0; 
  T     = (Mjd_UT1-MJD_J2000)/36525.0; 

  gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1
          + (0.093104-6.2e-6*T)*T*T; // [s]

  return  pi2*Frac(gmst/Secs);       // [rad], 0..2pi

}