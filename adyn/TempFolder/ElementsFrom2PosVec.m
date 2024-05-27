//------------------------------------------------------------------------------
//
// Elements 
//
// Purpose:
//
//   Computes orbital elements from two given position vectors and 
//   associated times 
//
// Input/Output:
//
//   GM        Gravitational coefficient
//             (gravitational constant * mass of central body)
//   Mjd_a     Time t_a (Modified Julian Date)
//   Mjd_b     Time t_b (Modified Julian Date)
//   r_a       Position vector at time t_a
//   r_b       Position vector at time t_b
//   <return>  Keplerian elements (a,e,i,Omega,omega,M)
//               a      Semimajor axis 
//               e      Eccentricity 
//               i      Inclination [rad]
//               Omega  Longitude of the ascending node [rad]
//               omega  Argument of pericenter  [rad]
//               M      Mean anomaly  [rad]
//             at time t_a 
//
// Notes:
//
//   The function cannot be used with state vectors describing a circular
//   or non-inclined orbit.
//
//------------------------------------------------------------------------------

Vector Elements ( double GM, double Mjd_a, double Mjd_b, 
                  const Vector& r_a, const Vector& r_b )
{
  
  // Variables
  
  double  tau, eta, p;
  double  n, nu, E, u;
  double  s_a, s_b, s_0, fac, sinhH;
  double  cos_dnu, sin_dnu, ecos_nu, esin_nu;
  double  a, e, i, Omega, omega, M;
  Vector  e_a(3), r_0(3), e_0(3), W(3);

  // Calculate vector r_0 (fraction of r_b perpendicular to r_a) 
  // and the magnitudes of r_a,r_b and r_0

  s_a = Norm(r_a);  e_a = r_a/s_a;
  s_b = Norm(r_b); 
  fac = Dot(r_b,e_a); r_0 = r_b-fac*e_a;
  s_0 = Norm(r_0);  e_0 = r_0/s_0;
  
  // Inclination and ascending node 

  W     = Cross(e_a,e_0);
  Omega = atan2 ( W(0), -W(1) );                     // Long. ascend. node 
  Omega = Modulo(Omega,pi2);
  i     = atan2 ( sqrt(W(0)*W(0)+W(1)*W(1)), W(2) ); // Inclination        
  if (i==0.0) 
    u = atan2 ( r_a(1), r_a(0) );
  else 
    u = atan2 ( +e_a(2) , -e_a(0)*W(1)+e_a(1)*W(0) );
  
  // Semilatus rectum
  
  tau = sqrt(GM) * 86400.0*fabs(Mjd_b-Mjd_a);   
  eta = FindEta ( r_a, r_b, tau );
  p   = pow ( s_a*s_0*eta/tau, 2 );   

  // Eccentricity, true anomaly and argument of perihelion

  cos_dnu = fac / s_b;    
  sin_dnu = s_0 / s_b;

  ecos_nu = p / s_a - 1.0;  
  esin_nu = ( ecos_nu * cos_dnu - (p/s_b-1.0) ) / sin_dnu;

  e  = sqrt ( ecos_nu*ecos_nu + esin_nu*esin_nu );
  nu = atan2(esin_nu,ecos_nu);

  omega = Modulo(u-nu,pi2);

  // Perihelion distance, semimajor axis and mean motion
  
  a = p/(1.0-e*e);
  n = sqrt ( GM / fabs(a*a*a) );

  // Mean anomaly and time of perihelion passage

  if (e<1.0) {
    E = atan2 ( sqrt((1.0-e)*(1.0+e)) * esin_nu,  ecos_nu + e*e );
    M = Modulo ( E - e*sin(E), pi2 );
  }
  else 
  {
    sinhH = sqrt((e-1.0)*(e+1.0)) * esin_nu / ( e + e * ecos_nu );
    M = e * sinhH - log ( sinhH + sqrt(1.0+sinhH*sinhH) );
  }

  // Keplerian elements vector

  return Vector(a,e,i,Omega,omega,M);

}