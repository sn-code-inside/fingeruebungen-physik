% -------------------------------------------------------------------------
% AufprallgeschwindigkeitKomet.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Aufprallgeschwindigkeit von Kometen/Meteoriten auf
% der Erde als Funktion verschiedener Parameter
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
RE = 6378;              % Erdradius in km
RJ = 71492;             % Jupiterradius in km
AE = 149597870;         % AE in km
aJ = 5.204*AE;          % Halbachse Jupiter
G  = 6.67*10^(-20);     % Gravitationskonstante in km^3/kg/s^2
MS = 1.9884*10^(30);    % Sonnemasse in kg
ME = 5.9723*10^(24);    % Erdmasse in kg
MJ = 1.899*10^(27);     % Jupitermasse in kg

V2 = sqrt(2*G*ME/RE);   % 2. Kosmische Geschwindigkeit
V3 = sqrt(G*MS/AE);     % 3. Kosmische Geschwindigkeit
V2J = sqrt(2*G*MS/aJ);  % 2. Kosmische Geschwindigkeit Jupiter

eSL9 = 0.215;           % Exzentrizität ShoemakerLevy 9    
aSL9 = 6.88648*AE;      % Halbachse ShoemakerLevy 9 
iSL9 = 6;               % Bahnneigung in Grad

pE = 3+2*ME*AE/MS/RE;


header1='Aufprallgeschwindigkeit als Funktion Halbachse'; % Gl. (3.238)
a=linspace(0.1,4,1000);
for m=1:2
    ex(m) = 0.75+(m-1)*0.15;
    for n=1:4
        i(n)  = 60*(n-1);
        vA((m-1)*4+n,:)  = real(V3*sqrt(pE-(1./a+2*cosd(i(n)).*sqrt(a*(1-ex(m)^2)))));
        lgdtitle((m-1)*4+n,:) = sprintf(' e = %3.2f  i = %3d° ',ex(m), i(n));
    end
end
figure('Name',header1);
for m=1:2
    for n= 1:4
        hold on; 
        if m == 1 
            LStyle =':';
        else
            LStyle ='-';
        end
        plot(a,vA((m-1)*4+n,:),'Color',Colors((m-1)*4+n+1,:),'LineWidth',2,'LineStyle', LStyle);
    end
end
xlim([0 2]);
ylim([0 80]);
grid on;
xlabel('Halbachse a in AE')
ylabel('v_A in km/s');
legend(lgdtitle,'location','northwest','NumColumns',2);
legend boxoff;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 18, 'linewidth', 1);


header1='Aufprallgeschwindigkeit als Funktion Perihelabstand'; % Gl. (3.239)
q=linspace(0.05,2,1000);
for n=1:4
    i(n)  = 60*(n-1);
    vB(n,:)  = real(V3*sqrt(pE-2*cosd(i(n)).*sqrt(2*q)));
    lgdtitleB(n,:) = sprintf('i = %3d° ', i(n));
 end
figure('Name',header1);
for n= 1:4
    hold on; 
    plot(q,vB(n,:),'Color',Colors(n+1,:),'LineWidth',2);
end
xlim([0.5 1.5]);
ylim([0 80]);
grid on;
xlabel('Perihelabstand q in AE')
ylabel('v_A in km/s');
legend(lgdtitleB,'location','northwest','NumColumns',1);
legend boxoff;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 18, 'linewidth', 1);

% Berechnung Jupiter Gl. (3.237)

veJ = sqrt(2*G*MJ/RJ),
U2 = (G*MS/aJ)*(3-(aJ/aSL9+2*cosd(iSL9)*sqrt(aSL9*(1-eSL9*eSL9)/aJ)));
VAJ =sqrt(veJ^2+U2);
fprintf('\n Escape GeschwindigkeitJupiter: %6.4f km/s \n',U2);
fprintf('\n Term U2: %6.4f km/s \n',veJ);
fprintf('\n Aufprallgeschwindigkeit SL9 auf Jupiter: %6.4f km/s \n',VAJ);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------