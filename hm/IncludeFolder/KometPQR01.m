% -------------------------------------------------------------------------
% KometPQR01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die heliozentrisch-ekliptikalen Koordinaten der Kometen 
% Gaussche Vektoren
% Eingabe:
%   T        Zeit in Julianischem Datum 
%   T0n      Periheldurchgang in Julianischem Datum 
%   qPn      Periheldistanz in AE
%   ePn      Exzentrizität
%   GM       Produkt aus Gravitationskonstante und Sonnenmasse AE^3*d°-2
% Ausgabe:
%  v        Geschwindigkeit des Kometen in kartesischen Koordinaten
%  xout_ell xyz Koordinaten des Kometen in Ellipsenberechnung
%  vout_ell Geschwindigkeit des Kometen in Ellipsenberechnung
%  xout_par xyz Koordinaten des Kometen in Parabelnäherungsberechnung
%  vout_par Geschwindigkeit des Kometen in Parabelnäherungsberechnung
%  time     Zeitstützstellen in JD
%  name     Name des Kometen
% -------------------------------------------------------------------------

% Kepler-Orbits PlanetPQR:
function KometData=KometPQR01(GM, T, BaPa,k)
%--------------------------------------------------------------------------
  pi2=2*pi;
  ePn=BaPa.eP(k);
  qPn=BaPa.qP(k);
  T0n=BaPa.T0(k);
  Komet(k).Name = BaPa.Name(k);
  Komet(k).Time = T;
  Tn=(T-2451545.0);
  T0n=(T0n-2451545.0);
  delta=1-ePn;
  invax=delta/qPn;
  tau=sqrt(GM)*(Tn-T0n);
  % Mittlere Anomalie 
  M  = tau*sqrt(invax*invax*invax);
  Komet(k).M  = M;
    aPn=1/invax;
    kf=sqrt(GM/aPn);
    Ecc=EAnom(M,ePn);
    cosE=cos(Ecc);
    sinE=sin(Ecc);
    fac=sqrt((1-ePn)*(1+ePn));
    rho=1-ePn.*cosE;
    x=aPn*(cosE-ePn); 
    y=aPn*fac.*sinE; 
    nu=2*atan(sqrt((1+ePn)/(1-ePn)).*tan(Ecc/2));
    R=aPn.*rho;
    vx=-kf.*sinE./rho;
    vy=kf*cosE*fac./rho;
    Komet(k).xout_ell   = [x;y;R*0]; 
    Komet(k).vout_ell   = [vx; vy; 0.0*vx];
    Komet(k).xyz = [0;0;0]; 
    Komet(k).v   = [0;0;0];
    fac =0.5*ePn;
    E2=100;
    E20=0;
    iter=0;
    kf=sqrt(GM/(qPn*(1+ePn)));
    while abs(E2-E20)>eps
      iter=iter+1;
      E20 = E2;
      A = 1.5*sqrt(fac/(qPn*qPn*qPn)).*tau;  
      B = (sqrt(A.*A+1.0)+A).^(1/3);
      u  = B - 1.0./B;  
      u2 = u.*u;  
      E2 = u2*(1.0-ePn)./fac;
      [c1,c2,c3]=Stumpff(E2); 
      fac = 3.0*ePn.*c3;
    end
    R  = qPn*(1 + u2.*c2*ePn./fac );
    x=qPn*(1.0-u2.*c2./fac);
    y=qPn*sqrt((1.0+ePn)./fac).*u.*c1;
    z=0*x;
    vx=-kf*y./R;
    vy=kf*(x./R+ePn);
    Komet(k).xout_par =[x;y;z];
    Komet(k).vout_par =[vx; vy; 0.0*vx];
    KometData=Komet(k);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
