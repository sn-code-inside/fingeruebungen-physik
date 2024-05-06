% -------------------------------------------------------------------------
% Analemma02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Modifizierte und ergänzte (2020) Public Domain Software von                                                                 
% Octave (Matlab) Peter Stallinga, February 2012                
% -------------------------------------------------------------------------
% Analemmaberechnung (Näherungsrechnung) für verschiedene geographische 
% Positionen in Alt-Azimut-Darstellung für 24 h. 
% Refraktion und Tagbogenverlängerung vernachlässigt. 
% Es wird eine einfache auf der Mittelpunktsgleichung basierende 
% Keplerlösung für die Erdbahn benutzt. 
% (Äquinoktium und Ekliptik J2000)
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
maanden = [31,28,31,30,31,30,31,31,30,31,30,31]; % months
y = 2000;
% Position of: Oberkochen Germany
lati  = 48.8;    % Geograph. Breite
longi = 10.1;    % Geograph. Länge
% Position of: Greenwich UK
% lati  = 45;    % Geograph. Breite
% longi = 0;    % Geograph. Länge
% Position of: Singapore UK
% lati  = 1;    % Geograph. Breite
% longi = 104;    % Geograph. Länge
% Position of: Tromsö Norway
% lati  = 69.6;    % Geograph. Breite
% longi = 18.9;    % Geograph. Länge
midpoint = 0;
for i=1:7
    midpoint = midpoint + maanden(i);
end
midpoint = midpoint+29;

latlong = sprintf(' Breite %4.1f° Länge %4.1f° (Uhrzeit in UT)', lati, longi);

% Prt 1: 
% Berechnung des Analemmas:

for h=1:24 %Stunden des Tages
    for n=1:365  %Tage des Jahres
        d = n;
        m = 1;
        while (d>maanden(m))  % Fortlaufende Tage des Jahres 
                              % umgerechnet in Tage des Monats
            d = d - maanden(m); % Monatstag
            m = m+1;            % Monat
        end
        %Berechnet Alt/Azimuth
        [az, alt] = azialti(y, m, d, h, 0, lati, longi); 
        az0h(h, n) = az;
        alt0h(h, n)= alt;
    end
end

% Teil 2
% Berechnung der Schnittlinien des 21. eines jeden Monats mit dem Analemma

for m=1:12
    n=0;
    d = 21;
    h=25.00;
    mins=0;
    alt = 100;
    while (h>0)
        [az,alt] = azialti(y, m, d, h, mins, lati, longi);
%        if (alt>0)
            n=n+1;
            az21m(m,n)=az;
            alt21m(m,n)=alt;
%         end
        mins=mins-1;
        if (mins<0)
            mins=mins+60;
            h = h-1;
        end
    end
    bserieslen(m)=n;
end

% Graphische Ausgabe
% Ausgabe der Aanalemma
figure1=figure('Name','Analemma');

for h=1:24
      plot(az0h(h,:), alt0h(h,:),'Color',Colors(3,:),'LineWidth',2);
      if alt0h(h,midpoint) > 0 
          text(az0h(h,midpoint), alt0h(h,midpoint),sprintf(' %02d:00',h));
      end
      hold on;
end
% Ausgabe der Monatsschnittlinie mit dem Anlemma
for m=1:12
    if bserieslen(m)>0
      plot(az21m(m,1:bserieslen(m)), alt21m(m,1:bserieslen(m)),'Color',Colors(m,:));
      hold on;
    end
end
% Beschriftung der Monatsschnittlinie mit dem Anlemma
for j=1:2
 for m=1:6
    k =(j-1)*6;
    dt = datetime(2000,m+k,21);
    MaxAlt = max(alt21m(m+k,:))-(-1)^j;
    MyLabel= string(dt,'dd.MM.');
    if MaxAlt>0 text(0, MaxAlt,MyLabel ,'Color',Colors(m+k,:)); end
 end
end
axis([-180 180 0 70])
text(-180,-2,'OST', 'FontSize',18);
text( 175,-2,'WEST', 'FontSize',18);
xlabel('Azimut (°)');
ylabel('Höhe (°)');
title(strcat('Analemma bei Position : ',latlong));
grid on;
grid minor;
set (gca,'FontSize',16);
%--------------------------------------------------------------------------
% Funktionen


% Berechnung von Höhe und Azimuth
function [azi,alti] = azialti(y, m, d, h, mins, lati, longi) 
    h = h+mins/60.0;

    % Part I: Inspired by http://www.stargazing.net/kepler/sun.html
    % 1. Find the days before J2000.0 (d) 
       d =  367 * y - floor(7 * (y +floor( (m + 9)/12))/ 4) + floor(275 * m / 9) + d - 730531.5 + h / 24;
   
    % 2. Find the Mean Longitude (L) of the Sun (siehe Formel 3.65 im Buch) 
       L = wrapTo360(280.461 + 0.9856474 * d);
      
    % 3. Find the Mean anomaly (g) of the Sun
       g = 357.528 + 0.9856003 * d;
      
    % 4. Find the ecliptic longitude (lambda) of the sun
        % (Mittelpunktsgleichung)
       lambda = L + 1.915 * sind(g) + 0.020 * sind(2*g);
   
    % 5. Find the obliquity of the ecliptic plane (epsilon)
       epsilon = 23.439 - 0.0000004 * d;
      
    % 6. Find the Right Ascension (alpha) and Declination (delta) of Sun
        % Formel 3.4 im Buch
       Y = cosd(epsilon) * sind(lambda);
       X = cosd(lambda);
       %Rektaszension    
       alpha = rad2deg(atan2(Y, X));     
       %Deklintaion
       delta = asind(sind(epsilon)*sind(lambda));
      
   % Part II: From http://www.geoastro.de/elevaz/basics/index.htm
   % compute Sidereal time at Greenwich:
   % (Naherungsformel)
       T = d/36525.0;  % Zeit in Julianischen Jahrhunderten
       theta0 = wrapTo360(280.46061837 + 360.98564736629*d + 0.000387933*T*T - T*T*T/38710000.0);
       theta = theta0 + longi;
       tau = theta - alpha; %Stundenwinkel)
       beta = lati;
       % Formel 3.4 im Buch
       alt = asin(sind(beta)*sind(delta) + cosd(beta)*cosd(delta)*cosd(tau));
       az  =  atan2(- sind(tau) , (cosd(beta)*tand(delta) - sind(beta)*cosd(tau)));
       alti= rad2deg(alt);
       azi = wrapTo180(rad2deg(az)-180);
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
