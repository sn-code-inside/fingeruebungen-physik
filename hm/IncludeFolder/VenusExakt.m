% -------------------------------------------------------------------------
% VenusExakt.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Das vorliegende Programmaterial basiert in Teilen 
% auf C++ Programmstrukturen und der Beschreibung/Erläuterung in  
%
% "Astronomie mit dem Personal Computer"
%
% von Oliver Montenbruck und Thomas Pfleger (Springer 1999). 
% Genehmigung des Verlages und der Autoren liegt vor.
% -------------------------------------------------------------------------
% Berechnet die heliozentrischen Koordinaten des Venus.
% Alle Werte beziehen sich auf  Äquinoktium des Datums
%  Eingabe:
%    t          Zeit in Julianischen Jahrhunderten seit JD2000
%    VenData    Eingelesene Tabelle Venus
%  Ausgabe:
%   VenPos.ekl  Heliozentrisch-ekliptikale Koordinaten Venus
%   VenPos.xyz  xyz Koordinaten (heliozentr.) Venus 
% -------------------------------------------------------------------------

function VenPos = VenusExakt(t,VenData)

pi2     = 2*pi;
Arcs    = 3600*180/pi;
dl = 0; dr = 0; db = 0;

% Mittlere Anomalien
M1 = pi2 * Frac ( 0.4861431 + 415.2018375*t );
M2 = pi2 * Frac ( 0.1400197 + 162.5494552*t );
M3 = pi2 * Frac ( 0.9944153 +  99.9982208*t );
M4 = pi2 * Frac ( 0.0556297 +  53.1674631*t );
M5 = pi2 * Frac ( 0.0567028 +   8.4305083*t );
M6 = pi2 * Frac ( 0.8830539 +   3.3947206*t );

% Störungen durch Merkur
for k =1:4
  iT = VenData(k,4);
  i5 = VenData(k,2);
  i6 = VenData(k,3);
  ct  = cos(i5*M2+i6*M1); st=sin(i5*M2+i6*M1);
  dlt = VenData(k,5)*ct+VenData(k,6)*st;
  drt = VenData(k,7)*ct+VenData(k,8)*st;
  dbt = VenData(k,9)*ct+VenData(k,10)*st;
  dl  = dl + dlt;
  dr  = dr + drt;
  db  = db + dbt; 
end
% Störungen durch Erde
for k =6:32
  iT = VenData(k,4);
  i5 = VenData(k,2);
  i6 = VenData(k,3);
  ct  = cos(i5*M2+i6*M3); st=sin(i5*M2+i6*M3);
  if ~(iT == 0) 
    ct =ct.*t.^iT; st = st.*t.^iT;  
  end
  dlt = VenData(k,5)*ct+VenData(k,6)*st;
  drt = VenData(k,7)*ct+VenData(k,8)*st;
  dbt = VenData(k,9)*ct+VenData(k,10)*st;
  dl  = dl + dlt;
  dr  = dr + drt;
  db  = db + dbt; 
end
% Störungen durch Mars
for k =34:36
  i5 = VenData(k,2);
  i7 = VenData(k,3);
  ct  = cos(i5*M2+i6*M4); st=sin(i5*M2+i7*M4);
  dlt = VenData(k,5)*ct+VenData(k,6)*st;
  drt = VenData(k,7)*ct+VenData(k,8)*st;
  dbt = VenData(k,9)*ct+VenData(k,10)*st;
  dl  = dl + dlt;
  dr  = dr + drt;
  db  = db + dbt; 
end
% Störungen durch Jupiter
for k =38:43
  i5 = VenData(k,2);
  i7 = VenData(k,3);
  ct  = cos(i5*M5+i6*M2); st=sin(i5*M5+i7*M2);
  dlt = VenData(k,5)*ct+VenData(k,6)*st;
  drt = VenData(k,7)*ct+VenData(k,8)*st;
  dbt = VenData(k,9)*ct+VenData(k,10)*st;
  dl  = dl + dlt;
  dr  = dr + drt;
  db  = db + dbt; 
end
% Störungen durch Saturn
for k =45:46
  i5 = VenData(k,2);
  i7 = VenData(k,3);
  ct  = cos(i5*M6+i6*M2); st=sin(i5*M6+i7*M2);
  dlt = VenData(k,5)*ct+VenData(k,6)*st;
  drt = VenData(k,7)*ct+VenData(k,8)*st;
  dbt = VenData(k,9)*ct+VenData(k,10)*st;
  dl  = dl + dlt;
  dr  = dr + drt;
  db  = db + dbt; 
end

% Ekliptikale Koordinaten ([rad],[AE])
dl = dl + 2.74*sin(pi2*(0.0764+0.4174*t))+0.27*sin(pi2*(0.9201+0.3307*t));
dl = dl + 1.9 + 1.8*t;
  
l = pi2 * Frac (0.3654783 + M2/pi2 + ( (5071.2+1.1.*t).*t + dl ) / 1296.0e3);
r = 0.7233482 - 0.0000002 * t +  dr * 1.0e-6;
b = ( -67.70 + ( 0.04 + 0.01.*t).* t  +  db ) / Arcs;

VenPos.ekl = [r; l; b];
VenPos.xyz = CalcXYZfromAngles(VenPos.ekl); % xyz helio-ekliptikal
VenPos.Name ="Venus";
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------