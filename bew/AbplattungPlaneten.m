% -------------------------------------------------------------------------
% AbplattungPlaneten.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Abplattung der Erde und von Gasplaneten
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

%%
% Anfangsbedingungen, Parameter
RE = 6378000;                       % mittlerer Erddurchmesser in m
ME = 5.974e24;                      % Masse Erde in kg
g0 = 9.81;                          % Schwerebeschleunigung
G  = 6.67384e-11;                   % Gravitationskonstante 
g0 = G*ME/RE/RE;
% Kreisrotationsfrequenz Erde
omegaE0 = 2*pi/86164;  
maxz = 1000;
% Geographische Breite
phi = linspace(0,90,maxz);
theta = 90-phi;


%%
% Aufgabenteil 1
fac = 6;
omegaE = 1*omegaE0;
geff = g0-RE*omegaE*omegaE.*cosd(phi).*cosd(phi);
gstr = g0*sqrt(1+RE^2*omegaE^4.*cosd(phi).*cosd(phi)*(1-2*g0/RE/omegaE^2)/g0^2);
figure();
plot(theta,geff,'LineStyle','-','LineWidth',1,'Color',Colors(2,:));
hold on;
plot(theta,gstr,'LineStyle',':','LineWidth',1,'Color',Colors(2,:));
omegaE = fac*omegaE0;
geff1 = g0-RE*omegaE*omegaE.*cosd(phi).*cosd(phi);
gstr1 = g0*sqrt(1+RE^2*omegaE^4.*cosd(phi).*cosd(phi)*(1-2*g0/RE/omegaE^2)/g0^2);
plot(theta,geff1,'LineStyle','-','LineWidth',1,'Color',Colors(3,:));
plot(theta,gstr1,'LineStyle',':','LineWidth',1,'Color',Colors(3,:));
lgd=legend ("g_{eff}","g'","g_{eff} 6\omega","g'    6\omega",'location', 'southwest');
set(lgd,'FontSize',12)
legend box off;
grid on;
xlabel('Winkel vom Nordpol in °','FontSize',14);
ylabel('g in m/s²','FontSize',14);
xlim([0,90]);
set(gca,'Fontsize', 16);

omegaE = 1*omegaE0;
alpha = atan2d(RE*omegaE*omegaE*sind(phi).*cosd(phi),g0-RE*omegaE*omegaE*cosd(phi).*cosd(phi));
figure();
plot(theta,alpha,'LineStyle','-','LineWidth',1,'Color',Colors(2,:));
hold on;
omegaE = fac*omegaE0;
alpha = atan2d(RE*omegaE*omegaE*sind(phi).*cosd(phi),g0-RE*omegaE*omegaE*cosd(phi).*cosd(phi));
plot(theta,alpha,'LineStyle','-','LineWidth',1,'Color',Colors(3,:));
lgd=legend ("\alpha","\alpha 6\omega",'location', 'northeast');
set(lgd,'FontSize',14)
legend box off;
grid on;
xlabel('Winkel vom Nordpol in °','FontSize',14);
ylabel("Winkelabwichung von g' und g_{eff}",'FontSize',14);
xlim([0,90]);
set(gca,'Fontsize', 16);

%%
% Aufgabenteil 2

omegaE = 1*omegaE0;

D   = RE^2*omegaE^2/6/g0;
h   = D*(1-3*(sind(phi)).^2);
figure();
plot(phi,h/1000);
plot(theta,h/1000,'LineStyle','-','LineWidth',2,'Color',Colors(3,:));

RAex = 6378388;
RPex = 6356752;
RA = RE + h(1);
RP = RE + h(maxz);
DR = (RA - RP)/1000;
xel = linspace(RA,0,1000);
yel = sqrt((1-xel.*xel/RA^2))*RP;
rel = sqrt(xel.^2+yel.^2);
rhoE = ME/(4*pi*RE^3/3);
r = (RE + h);

hold on;
plot(theta,(rel-r)/1000,'LineStyle',':','LineWidth',2,'Color',Colors(3,:));
lgd=legend ("von mittlerer Sphäre","von Rotationsellipsoid",'location', 'southeast');
set(lgd,'FontSize',14)
legend box off;
grid on;
xlabel('Winkel vom Nordpol in °','FontSize',14);
ylabel("Abweichung in km ",'FontSize',14);
xlim([0,90]);
set(gca,'Fontsize', 16);

%%
% Aufgabenteil 3

kappa1_E = 1+omegaE^2*RA^3/2/G/ME;
% Erde
func=@(kap)(kap^2+2)*(kap^2-1)^(-3/2)*atan(sqrt(kap^2-1))-3/(kap^2-1)- ...
    omegaE^2/2/pi/G/rhoE;
kappa2_E = fzero(func,[0.5 1.5]);
RP1 = kappa1_E*RA;
RP2 = kappa2_E*RA;
DR1 = RA - RP1;
DR2 = RA - RP2;
gamma = DR2/DR1;
fprintf("\n Berechnung der Differenz Äquator - Polradius \n");
fprintf("\n Planet  \t kappa1  \t kappa2 \t    gamma  \t Differenz1\t Differenz2 \tExperiment     ");
fprintf("\n");
fprintf('\n Erde    \t %8.6f \t %8.6f \t %8.3f \t %4.0f km \t  %4.0f km \t \t %4.0f km', ...
    kappa1_E, kappa2_E, gamma, -DR1/1000, -DR2/1000, -(RPex-RAex)/1000);
fprintf("\n");

omegaP = [1.74533E-04;	1.63115E-04;	1.02666E-04;	1.08406E-04];
rhoP   = [1326;	687;	1271;	1638];
RAP    = [71492000;	60268000;	25559000;	24764000];
MP     = [1.89810E+27;	5.68296E+26;	8.68100E+25;	1.02403E+26];
DPPex  = [133708;	108728;	49946;	48682];
RPPex  = DPPex*1000/2;

kappa1_P = 1+omegaP.^2.*RAP.^3/2/G./MP;

myfunc=@(kap, c1, c2)(kap^2+2)*(kap^2-1)^(-3/2)*atan(sqrt(kap^2-1))-3/...
    (kap^2-1)-c1^2/2/pi/G/c2;
for k=1:4
    c1 = omegaP(k);
    c2 = rhoP(k);
    func = @(x) myfunc(x,c1, c2);    % function of x alone
    kappa2_P(k) = fzero(func,[0.5 1.5]);
end
G = 6.67384E-11;

RPP1 = kappa1_P.*RAP;
RPP2 = kappa2_P.*RAP;
DRP1 = RAP - RPP1;
DRP2 = RAP - RPP2;
gamma = DRP2./DRP1;
Planet =["Jupiter"; "Saturn "; "Uranus "; 'Neptun '];
for k=1:4
fprintf('\n %s \t %8.6f \t %8.6f \t %8.3f \t %4.0f km \t  %4.0f km \t \t %4.0f km', ...
    Planet(k,:), kappa1_P(k), kappa2_P(k), gamma(k), -DRP1(k)/1000, ...
    -DRP2(k)/1000, -(RPPex(k)-RAP(k))/1000);
end

fprintf('\n');
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
