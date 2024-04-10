% -------------------------------------------------------------------------
% FindRiseSet.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm entwicklt nach Code in Montenbruck
% O. Montenbruck and T. Pfleger
% Astronomie mit dem Personal Computer
% Springer, Berlin 1999
% -------------------------------------------------------------------------
% Nullstellensuche Sonnen-Auf- und Untergang im Intervall 0..24 h mit
% quadratischer Interpolation unter Berücksichtung von Polarnacht und -tag
% Verwendung in SunRise02.m 
% -------------------------------------------------------------------------

function [xTRise, xTSet, xRise, xSet, xAbove] = FindRiseSet(MJT, lambda, phi)
% Eingabe: 
%   MJT     =  Modifiziertes julianisches Datum 
%   lambda  =  geogr. Länge     [lambda]    =  °
%   phi     =  geogr. Breite    [phi]       = rad 
% Ausgabe: 
%   xTRise  =  Zeitpunkt Sonnenaufgang
%   xTSet   =  Zeitpunkt Sonnenuntergang
%   xRise   =  Sonnenaufgang ja/nein
%   xSet    =  Sonnenuntergang ja/nein
%   xAbove  =  Polarnacht/"normaler" Tag/Polartag 

% Berücksichtigung der Refraktion zur Ermittlung der scheinbaren Höhe 
   sinhr  =sind(-50/60);
   
   sinphi =sin(phi);
   cosphi=cos(phi);
   band=24;         % Taglänge 
   stst= 48;        % Stützstellen 
   x_s=0;           % Startwert Stunde
   dx=band/stst/2;  % zeitliche Schrittweite 
   
   % Berechnung des Sinus der scheinbaren Sonnenhöhe 
   y_minus=SinAlt(MJT,x_s,lambda,cosphi,sinphi)-sinhr;
   if y_minus > 0 
       Above = 1;   % Sonne über dem Horizont  
   else 
       Above = 0;   % Sonne unter dem Horizont 
   end
   
   % Initialisierung 
   Rise = false;
   Set = false;
   TRise=0;
   TSet=0;
   
   % Iteration über die Tageszeit 
   for n=1:stst
    % Berechnung dreier benachbarter Zeitwerte
    x_minus = x_s+(2*n-2)*dx;
    x_0     = x_minus+dx;
    x_plus  = x_0+dx;
    % Berechnung des Sinus der scheinbaren Sonnenhöhe für alle drei
    % Zeitpunkte
    y_minus=SinAlt(MJT,x_minus,lambda,cosphi,sinphi)-sinhr;
    y_0=(SinAlt(MJT,x_0,lambda,cosphi,sinphi)-sinhr);
    y_plus=(SinAlt(MJT,x_plus,lambda,cosphi,sinphi)-sinhr);
    
    [xroot1, xroot2, xnroot]= squintpol(y_minus,y_0,y_plus);
    if xnroot == 1
      if y_minus < 0 
              TRise= x_0+dx*xroot1;
              Rise=true;
      else
              TSet= x_0+dx*xroot2;
              Set=true;
      end
    end
    if xnroot==2 
     if y_0 < 0 %Ueberpruefen, sonst ye aus Squintpol uebergeben lassen
              TSet  = x_0+dx*xroot1;
              TRise = x_0+dx*xroot2;
     else
              TSet  = x_0+dx*xroot2;
              TRise = x_0+dx*xroot1;
     end
     Set=true;
     Rise=true;
    end
   end
   xTRise = TRise;
   xTSet = TSet;
   xRise = Rise;
   xSet = Set;
   % für Polarregionen 
   if (xRise==false) && (xSet==false)
       % Polartag
       if y_minus > 0  
           xAbove = 1;
       % Polarnacht 
       else
           xAbove = -1;
       end
   else
       xAbove= 0;
   end   
end      

% -------------------------------------------------------------------------
% SinAlt
% Berechnung des Sinus der Höhe der Sonne
% Verwendung in FindRiseSet (siehe oben)  
% -------------------------------------------------------------------------

function sinusalt=SinAlt(MJT, hour, lambda, cosphi, sinphi)
% Eingabe: 
%   MJT     = modifiziertes julianisches Datum 
%   hour    = Zeitpunkt in Stunden 
%   lambda  = geogr. Länge 
%   cosphi  = Cosinus der geogr. Breite
%   sinphi  = Sinus der geogr. Breite
% Ausgabe: 
%   sinusalt = Sinus der Sonnenhöhe 

  MJD=0;
  MJD=MJT+hour/24;              % Umrechnung der Zeit 
  eps=deg2rad(23.43929111);     % Schiefe der Ekliptik 
  T=MJD+2400000.5;
  [xRA,xDec]=KeplerSonne(T,eps);  % Keplerlösung f. Rektaszension/Deklination 
  tau=GMST(MJD)+deg2rad(lambda)-xRA;
  tau = wrapToPi(tau);          % Berechnung des Stundenwinkels 
  sinusalt=sinphi*sin(xDec)+cosphi*cos(xDec)*cos(tau);
end

% -------------------------------------------------------------------------
% squintpol 
% Quadratische Interpolation
% Verwendung in FindRiseSet (siehe oben)  
% -------------------------------------------------------------------------

function [xroot1, xroot2, xnroot]= squintpol(y_minus,y_0,y_plus)
% Eingabe: 
%   y_minus, y_0, y_plus = Sinus der Höhe der Sonne für drei Zeitpunkte 
% Ausgabe: 
%   xroot1      = 1. Wurzel 
%   xroot2      = 2. Wurzel 
%   xnroot      = Zahl der Wurzeln 
   nroot = 0;
   xnroot = nroot;
   root1 = 0;
   root2 = 0;
   a = 0.5*(y_plus+y_minus)-y_0;
   b = 0.5*(y_plus-y_minus);
   c = y_0;
   rd = 1-4*a*c/(b*b);
   if rd ==1  % Es gibt genau eine Wurzel !
       root1=-b/(2*a);
       root2=-root1;
       nroot=1;
   else
     if rd >=0 % Es gibt eine oder zwei Wurzeln !
       root1=(-2*c/b)/(1-sqrt(rd));
       root2=(-2*c/b)/(1+sqrt(rd));
       if abs(root1)<=1 
           nroot=nroot+1;
       end
       if abs(root2)<=1 
          nroot=nroot+1;
       end
       if (root1 < -1.0) || (root1 > 1.0)
           root1=root2;
       end
     end
   end
   xroot1=root1;
   xroot2=root2;
   xnroot=nroot;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------