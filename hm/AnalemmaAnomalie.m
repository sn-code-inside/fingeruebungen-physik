% -------------------------------------------------------------------------
% AnalemmaAnomalie.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Dank an Rudolf von Buenau und Kenneth von Buenau fuer Unterstuetzung und
% wertvolle ANregungen.
%-------------------------------------------------------------------------
% Programm berechnet die Analemmas auf unseren Planeten direkt aus 
% der exzentrischen Anomalie. Man erhält die AnalemmaForm, aber da man die 
% Keplergleichung nicht loest, keine Zeitabhaengigkeit, also auch keine
% Zeitgleichung.
% 
% Für die Parameter der Planeten wird das Aequinoktium J2000 benutzt.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
%opengl('save', 'software');

% Initialisierung
% Exzentrizitaet
exzP = [0.2056,0.0068,0.0167,0.0934,0.0485,0.0555,0.0464,0.0095,0.2488];  
% Polneigung in °
epsP = [0.02,2.64,23.44,25.19,3.13,26.73,82.23,28.23,57.47];   
% Planetozentrische Laenge des Perihels in ° 
% Selbst berechnet oder von Michael Allison und Megan McEwan
% Planetary and Space Science 48(2000) 215  
perP = [8.0,253.8,282.94,251.0,57.1,279.5,185.5,2.2,184.5];   %eigene Ber.
perP = [325.7,253.8,282.94,251.0,57.1,279.5,185.5,2.2,184.5]; %Allison
Name = ["Merkur" "Venus" "Erde" "Mars" "Jupiter" "Saturn" "Uranus"...
        "Neptun" "Pluto"];

% Berechnung und direkte Ausgabe
figure()
for PlanetNr = 1:9
    [tau, dd] = analemma (exzP(PlanetNr), epsP(PlanetNr), perP(PlanetNr));
    subplot(3,3,PlanetNr);
    plot(tau,dd,'color',Colors(PlanetNr,:),'LineWidth',2);
    am=max(abs(tau)); % window size x
    am = round(1.2*am); 
    om=max(abs(dd));  % window size y
    om = round(1.2*om); 
    if om <1 
        om=0.05; 
    end
    axis([-am, am, -om, om])
    title(Name(PlanetNr));
    if PlanetNr > 6 
        xlabel('ZGL in min'); 
    end
    if PlanetNr==1 || PlanetNr==4 || PlanetNr==7 
        ylabel('Deklination \delta in °'); 
    end
    grid on;
    set(gca,'Fontsize',16);
 end

%%
% Berechnung des Analemma über die exzentrische Anomalie E

function [tau, d_WS] = analemma (exzp, epsp, perp)  
% Berechnung im Bahn-Koordinatensystem
    E = 0:360;                        % exzentrische Anomalie f. 1 Jahr
    M = E - exzp * sind(E) * 180/pi;  % mittlere Anomalie, Keplergleichung
    % wahre Anomalie aus Formel (46) Kap. Himmelsmechanik
    ups = 2*atan2d(sqrt((1+exzp)/(1-exzp)),1./tand(E/2)); 

% Übergang ins planetozentrische KOS
    l_WS = ups + perp;                % Laenge wahre Sonne (WS)
    l_MS = M + perp;                  % Laenge mittlere Sonne (MS)
 
% Übergang ins planetozentrische-äquatoriale KOS (analog GA(sf)-KOS)
    % Tabelle 7 Kap. Himmelmechanik
    alpha_WS = atand(cosd(epsp)*tand(l_WS));        % Rektaszension WS
    % Abbildung in einem über über 360° zusammenhängenden Intervall 
    alpha_WS = alpha_WS + floor((l_WS+90)/180)*180;
    alpha_MS = l_MS;                                % Rektaszension MS2        
    d_WS = asind(sind(epsp) * sind(l_WS));          % Deklination WS
% Stundenwinkel Wahre Sonne, Formel (81) Kap. Himmelsmechanik
    tau = (alpha_MS - alpha_WS) * 4;                % Stundenwinkel in min
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
