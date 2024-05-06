% -------------------------------------------------------------------------
% PlanetParameterEinlesen.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Einlesen der Orbitalparameter nach Montenbruck/Pfleger
% -------------------------------------------------------------------------

function [aP,epsP,iP,M0P,nP,OmegaP,wP,L0P,TP,PlanetTitle]=PlanetParameterEinlesen
PlanetTitle   =  strings([10,25]);
PlanetTitle(1,:)= ' Merkur';
PlanetTitle(2,:)= ' Venus';
PlanetTitle(3,:)= ' Erde';
PlanetTitle(4,:)= ' Mars';
PlanetTitle(5,:)= ' Jupiter';
PlanetTitle(6,:)= ' Saturn';
PlanetTitle(7,:)= ' Uranus';
PlanetTitle(8,:)= ' Neptun';
PlanetTitle(9,:)= ' Pluto';
PlanetTitle(10,:)= ' Sonne';
aP=zeros(10);
epsP=zeros(10);
iP=zeros(10);
OmegaP=zeros(10);
wP=zeros(10);
L0P=zeros(10);
M0P=zeros(10);
TP=zeros(10);

% case Mercury = 1: 
      aP(1) =  0.387099;
      epsP(1) = 0.205634; 
      M0P(1) = 174.7947; 
      nP(1) = 149472.6738;
      OmegaP(1) =  48.331;   
      iP(1) = 7.0048;   
      wP(1) =  77.4552; 
      TP(1) = 0.0;
% case Venus = 2: 
      aP(2) =  0.723332;
      epsP(2) = 0.006773; 
      M0P(2) =  50.4071; 
      nP(2) = 58517.8149;
      OmegaP(2) =  76.680;   
      iP(2) = 3.3946;  
      wP(2) = 131.5718;
      TP(2) = 0.0;
% case Earth = 3:   
      aP(3) =  1.000000; 
      epsP(3) = 0.016709; 
      M0P(3) = 357.5256; 
      nP(3) = 35999.3720;
      OmegaP(3) = 174.876;   
      iP(3) = 0.0000;   
      wP(3) = 102.9400; 
      TP(3) = 0.0;
% case Mars = 4:  
      aP(4) =  1.523692; epsP(4) = 0.093405; M0P(4) =  19.3879; nP(4) = 19140.3023;
      OmegaP(4) =  49.557;   iP(4) = 1.8496;   wP(4) = 336.0590; TP(4) = 0.0;
% case Jupiter = 5: 
      aP(5) =  5.204267; epsP(5) = 0.048775; M0P(5) =  18.8185; nP(5) =  3033.6272;
      OmegaP(5) = 100.4908;  iP(5) = 1.3046;   wP(5) =  15.5576; TP(5) = 0.0;
% case Saturn = 6:  
      aP(6) =  9.582018; epsP(6) = 0.055723; M0P(6) = 320.3477; nP(6) =  1213.8664;
      OmegaP(6) = 113.6427;  iP(6) = 2.4852;   wP(6) =  89.6567; TP(6) = 0.0;
% case Uranus = 7: 
      aP(7) = 19.229412; epsP(7) = 0.044406; M0P(7) = 142.9559; nP(7) =   426.9282;
      OmegaP(7) =  73.9893;  iP(7) = 0.7726;   wP(7) = 170.5310; TP(7) = 0.0;
% case Neptune = 8:
      aP(8) = 30.103658; epsP(8) = 0.011214; M0P(8) = 267.7649; nP(8) =   217.9599;
      OmegaP(8) = 131.7942;  iP(8) = 1.7680;   wP(8) =  37.4435; TP(8) = 0.0;
% case Pluto = 9:
      aP(9) = 39.264230; epsP(9) = 0.244672; M0P(9) =  15.0233; nP(9) =   146.3183;
      OmegaP(9) = 110.2867;  iP(9) = 17.1514;  wP(9) =  224.0499; TP(9) = 0.0;
% case Sun = 10:
      aP(10) = 0; epsP(10) = 1; 
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
