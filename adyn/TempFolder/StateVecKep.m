//------------------------------------------------------------------------------
//
// State
//
// Purpose:
//
//   Computes the satellite state vector from osculating Keplerian elements 
//   for elliptic orbits
//
// Input/Output:
//
//   GM        Gravitational coefficient
//             (gravitational constant * mass of central body)
//   Kep       Keplerian elements (a,e,i,Omega,omega,M) with
//               a      Semimajor axis 
//               e      Eccentricity 
//               i      Inclination [rad]
//               Omega  Longitude of the ascending node [rad]
//               omega  Argument of pericenter  [rad]
//               M      Mean anomaly at epoch [rad]
//   dt        Time since epoch
//   <return>  State vector (x,y,z,vx,vy,vz)
//
// Notes:
//
//   The semimajor axis a=Kep(0), dt and GM must be given in consistent units, 
//   e.g. [m], [s] and [m^3/s^2]. The resulting units of length and velocity  
//   are implied by the units of GM, e.g. [m] and [m/s].
//
//------------------------------------------------------------------------------

Vector State ( double GM, const Vector& Kep, double dt )
{

  // Variables

  double  a,e,i,Omega,omega,M,M0,n;
  double  E,cosE,sinE, fac, R,V;
  Vector  r(3),v(3);
  Matrix  PQW(3,3);

  // Keplerian elements at epoch
  
  a = Kep(0);  Omega = Kep(3);
  e = Kep(1);  omega = Kep(4); 
  i = Kep(2);  M0    = Kep(5);

  // Mean anomaly

  if (dt==0.0) {
    M = M0;
  }
  else {
    n = sqrt (GM/(a*a*a));
    M = M0 +n*dt;
  };

  // Eccentric anomaly
  
  E  = EccAnom(M,e);   

  cosE = cos(E); 
  sinE = sin(E);

  // Perifocal coordinates

  fac = sqrt ( (1.0-e)*(1.0+e) );  

  R = a*(1.0-e*cosE);  // Distance
  V = sqrt(GM*a)/R;    // Velocity

  r = Vector ( a*(cosE-e), a*fac*sinE , 0.0 );
  v = Vector ( -V*sinE   , +V*fac*cosE, 0.0 ); 

  // Transformation to reference system (Gaussian vectors)
  
  PQW = R_z(-Omega) * R_x(-i) * R_z(-omega);

  r = PQW*r;
  v = PQW*v;

  // State vector 
  
  return Stack(r,v);

}
