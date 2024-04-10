% -------------------------------------------------------------------------
% JupiterExakt.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
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
% JupPos = JupiterExakt(t,JupData)
% Berechnet die heliozentrischen Koordinaten des Jupiter.
% Alle Werte beziehen sich auf  Äquinoktium des Datums
%  Eingabe:
%   t         	Zeit in Julianischen Jahrhunderten seit JD2000
%   JupData     Eingelesene Tabelle Jupiter
%  Ausgabe:
%   JupPos.ekl  Heliozentrisch-ekliptikale Koordinaten Jupiter
%   JupPos.xyz  xyz Koordinaten (heliozentr.) Jupiter 
% -------------------------------------------------------------------------

function JupPos = JupiterExakt(t,JupData)

pi2     = 2*pi;
Arcs    = 3600*180/pi;
dl = 0; dr = 0; db = 0;

% Mittlere Anomalien
M5 = pi2 * Frac ( 0.0565314 + 8.4302963*t );
M6 = pi2 * Frac ( 0.8829867 + 3.3947688*t );
M7 = pi2 * Frac ( 0.3969537 + 1.1902586*t );
% Störungen durch Saturn
for k =1:48
  iT = JupData(k,4);
  i5 = JupData(k,2);
  i6 = JupData(k,3);
  ct  = cos(i5*M5+i6*M6); st=sin(i5*M5+i6*M6);
  if ~(iT == 0) 
    ct =ct.*t.^iT; st = st.*t.^iT;  
  end
  dlt = JupData(k,5)*ct+JupData(k,6)*st;
  drt = JupData(k,7)*ct+JupData(k,8)*st;
  dbt = JupData(k,9)*ct+JupData(k,10)*st;
  dl  = dl + dlt;
  dr  = dr + drt;
  db  = db + dbt; 
end
% Störungen durch Uranus
for k =50:51
  i5 = JupData(k,2);
  i7 = JupData(k,3);
  ct  = cos(i5*M5+i6*M6); st=sin(i5*M5+i7*M7);
  dlt = JupData(k,5)*ct+JupData(k,6)*st;
  drt = JupData(k,7)*ct+JupData(k,8)*st;
  dbt = JupData(k,9)*ct+JupData(k,10)*st;
  dl  = dl + dlt;
  dr  = dr + drt;
  db  = db + dbt; 
end
% Störungen durch Saturn und Uranus
phi = (2*M5-6*M6+3*M7); ct=cos(phi); st=sin(phi);
dl = dl - 0.8*ct + 8.5*st;
dr = dr - 0.1*ct;

phi = (3*M5-6*M6+3*M7); ct=cos(phi); st=sin(phi);
dl  = dl + 0.4*ct + 0.5*st;
dr  = dr - 0.7*ct + 0.5*st;
db  = db - 0.1*ct;

% Ekliptikale Koordinaten ([rad],[AE])
l = pi2 * Frac (0.0388910 + M5/pi2 + ( (5025.2+0.8.*t).*t + dl ) / 1296.0e3);
r = 5.208873 + 0.000041.*t  +  dr * 1.0e-5;
b = ( 227.3 - 0.3.*t + db ) / Arcs;
JupPos.ekl = [r; l; b];
JupPos.xyz = CalcXYZfromAngles(JupPos.ekl); % xyz helio-ekliptikal
JupPos.Name = "Jupiter";
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------