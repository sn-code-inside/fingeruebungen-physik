% -------------------------------------------------------------------------
% SatPQR_perturbed_orbit.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Astrodynamik aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die geozentrisch-äquatorialen Koordinaten von Satelliten aus
% den Bahnpaarmetern fÃ¼r eine bestimmte Zeit unter Berücksichtigung
% säkularer Störungen
%
% Eingabe:
%   T                               = Zeit in Julianischem Datum
%   aPn,ePn,iPn,MPn,OmegaPn,omegaPn = Orbitparameter in km
%                                     bzw. rad Ã¼bergeben
%   OmegaDot, omegaDot, nDot        = säkulare Störungen in rad/day
% Ausgabe:
%  xyz      = geozentrisch-äquatorialen Koordinaten (xyz) des Satelliten
%  vel      = Gesschwindigkeitskomponenten (xyz) des Satelliten
%  Sat.az; Sat.el, Sat.r   = geozentrisch-äquatorialen Koordinaten
%                            (alpha, delta, r) des Satelliten
%  Time     = ZeitstÃ¼tzstellen in JD
%  Name     = Name des Satelliten
% -------------------------------------------------------------------------

function SatData= SatPQR_perturbed_orbit(T, SatSet, GME, OmegaDotd, omegaDotd, nDotd)
  % Zeit  T in Julianischen Datum/Tagen 
  
  aP = SatSet.BaPa(1);
  eP = SatSet.BaPa(2);
  iP = SatSet.BaPa(3);
  MP = SatSet.BaPa(4);  %Mittlere Anomalie zum Startzeitpunkt
  OmegaP0 = SatSet.BaPa(5);
  omegaP0 = SatSet.BaPa(6);
  nP  = sqrt(GME/aP^3)*86400;  %Umlaufgeschwindigkeit in 1/day
% Mittlere Anomalie
  % Berechnung der exzentrischen Anomalie mit Startwert M(1) = MPn
  for k=1:length(T)
      M(k) = MP + (nP + nDotd + omegaDotd)*(T(k)-T(1)) ;
      Ecc(k)=EAnom(M(k),eP);
      OmegaP(k) = OmegaP0+OmegaDotd*(T(k)-T(1));
      omegaP(k) = omegaP0+omegaDotd*(T(k)-T(1));
  end
% Berechnung der Bahndaten
  cosE=cos(Ecc);
  sinE=sin(Ecc);
  fac=sqrt((1-eP)*(1+eP));
  % Bahnkoordinaten 
  dis     = aP.*(1-eP*cosE);
  Vel     = sqrt(GME*aP)./dis;
  uvw     = [aP.*(cosE-eP);aP*fac.*sinE;0.*Ecc]; 
  vel     = Vel.*[-sinE;fac.*cosE;0.*Ecc];
  for k=1:length(T)
     M_PQR = PQR(-OmegaP(k), -iP,-omegaP(k));
     Sat.xyz(:,k) = mtimes(M_PQR,uvw(:,k));
     Sat.vel(:,k) = mtimes(M_PQR,vel(:,k));
  end
%   Sat.xyz = mtimes(M_PQR,uvw);
%   Sat.vel = mtimes(M_PQR,vel);
  xP = Sat.xyz(1,:);
  yP = Sat.xyz(2,:);
  zP = Sat.xyz(3,:);
 [Sat.az,Sat.el,Sat.r] = cart2sph(xP,yP,zP);
  Sat.Time = T; 
  Sat.Name = SatSet.Name;
  SatData=Sat;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------