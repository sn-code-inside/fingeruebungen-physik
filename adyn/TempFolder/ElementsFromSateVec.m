//------------------------------------------------------------------------------
//
// Elements
//
// Purpose:
//
//   Computes the osculating Keplerian elements from the satellite state vector
//   for elliptic orbits
//
// Input/Output:
//
//   GM        Gravitational coefficient
//             (gravitational constant * mass of central body)
//   y         State vector (x,y,z,vx,vy,vz)
//   <return>  Keplerian elements (a,e,i,Omega,omega,M) with
//               a      Semimajor axis 
//               e      Eccentricity 
//               i      Inclination [rad]
//               Omega  Longitude of the ascending node [rad]
//               omega  Argument of pericenter  [rad]
//               M      Mean anomaly  [rad]
//
// Notes:
//
//   The state vector and GM must be given in consistent units, 
//   e.g. [m], [m/s] and [m^3/s^2]. The resulting unit of the semimajor
//   axis is implied by the unity of y, e.g. [m].
//
//   The function cannot be used with state vectors describing a circular
//   or non-inclined orbit.
//
//------------------------------------------------------------------------------

Vector Elements ( double GM, const Vector& y )
{

  // Variables

  Vector  r(3),v(3),h(3);
  double  H, u, R;
  double  eCosE, eSinE, e2, E, nu;
  double  a,e,i,Omega,omega,M;

  r = y.slice(0,2);                                  // Position
  v = y.slice(3,5);                                  // Velocity
  
  h = Cross(r,v);                                    // Areal velocity
  H = Norm(h);

  Omega = atan2 ( h(0), -h(1) );                     // Long. ascend. node 
  Omega = Modulo(Omega,pi2);
  i     = atan2 ( sqrt(h(0)*h(0)+h(1)*h(1)), h(2) ); // Inclination        
  u     = atan2 ( r(2)*H, -r(0)*h(1)+r(1)*h(0) );    // Arg. of latitude   

  R  = Norm(r);                                      // Distance           

  a = 1.0 / (2.0/R-Dot(v,v)/GM);                     // Semi-major axis    

  eCosE = 1.0-R/a;                                   // e*cos(E)           
  eSinE = Dot(r,v)/sqrt(GM*a);                       // e*sin(E)           

  e2 = eCosE*eCosE +eSinE*eSinE;
  e  = sqrt(e2);                                     // Eccentricity 
  E  = atan2(eSinE,eCosE);                           // Eccentric anomaly  

  M  = Modulo(E-eSinE,pi2);                          // Mean anomaly

  nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          // True anomaly

  omega = Modulo(u-nu,pi2);                          // Arg. of perihelion 
 
  // Keplerian elements vector

  return Vector(a,e,i,Omega,omega,M);

}

