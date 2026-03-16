% -------------------------------------------------------------------------
% AbplattungPlaneten.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Abplattung der Erde und von Gasplaneten
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('IncludeFolder')
addpath('Data')
Colors = GetColorLines;

%%
% Anfangsbedingungen, Parameter
RAexE = 6378137.00 ;                 % Erd-Äquatorradius in m
RPexE = 6356752.00 ;                 % Erd-Polradius in m
RE    = 6371000;                     % mittlerer Erdradius in m
G     = 6.67384e-11;                 % Gravitationskonstante 
ME    = 5.974e24;                    % Masse Erde in kg
g0    = 9.81;                        % Schwerebeschleunigung
g0    = G*ME/RE/RE;
% Kreisrotationsfrequenz Erde
omegaE0 = 2*pi/86164;  
maxz    = 1000;
% Geographische Breite
phi     = linspace(0,90,maxz);
theta   = 90-phi;
kappaexp_E = RAexE/RPexE;
fexp_E     = 1/(1-1/kappaexp_E);

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

% Erde 

%Modell 1
fA1_E = (omegaE0^2*RA^3/2/G/ME)^(-1);
RP1 = (1-1/fA1_E)*RAexE;
DR1   = RAexE - RP1;

% Modell 2
fA2_E = (5*omegaE0^2*RE^3/4/G/ME)^(-1);
RP2   = (1-1/fA2_E)*RAexE;
DR2   = RAexE - RP2;

% Erde
func=@(kap)(kap^2+2)*(kap^2-1)^(-3/2)*atan(sqrt(kap^2-1))-3/(kap^2-1)- ...
    omegaE^2/2/pi/G/rhoE;
eta3_E  = fzero(func,[0.5 1.5]);
fA3_E   = 1/(1-1/eta3_E);
RP3     = (1-1/fA3_E)*RAexE;
DR3     = RAexE - RP3;
gamma = DR3/DR1;
fprintf("\n Berechnung Abplattung aus Differenz Äquator - Polradius \n");
fprintf("\n Planet    \t fA1     \t fA2      \t fA3     \t fA_exp \t gamma    \tDelta R_1 \tDelta R_2 \tDelta R_3 \tDelta R (real)    ");
fprintf("\n");
fprintf('\n Erde    \t %5.2f   \t %5.2f   \t %5.2f  \t %5.2f \t %5.3f  \t%4.1f km  \t%4.1f km \t%4.1f km \t%4.1f km', ...
    fA1_E, fA2_E, fA3_E, fexp_E, gamma,  DR1/1000, DR2/1000, DR3/1000, -(RPexE-RAexE)/1000);
fprintf("\n");


% Planeten
omegaP = [1.74533E-04;	1.63115E-04;	1.02666E-04;	1.08406E-04];
rhoP   = [1326;	687;	1271;	1638];
RAP    = [71492000;	60268000;	25559000;	24764000];
MP     = [1.89810E+27;	5.68296E+26;	8.68100E+25;	1.02403E+26];
DPPex  = [133708;	108728;	49946;	48682];
RPPex  = DPPex*1000/2;
fAexp_P= [15.41 10.21 43.62 58.54];

%Modell 1
fA1_P = (omegaP.^2.*RAP.^3/2/G./MP).^(-1);
RPP1  = (1-1./fA1_P).*RAP;
DRP1  = RAP - RPP1;

% Modell 2
fA2_P = (5*omegaP.^2.*RAP.^3/4/G./MP).^(-1);
RPP2  = (1-1./fA2_P).*RAP;
DRP2  = RAP - RPP2;

% Modell 3
myfunc=@(kap, c1, c2)(kap^2+2)*(kap^2-1)^(-3/2)*atan(sqrt(kap^2-1))-3/...
    (kap^2-1)-c1^2/2/pi/G/c2;
for k=1:4
    c1 = omegaP(k);
    c2 = rhoP(k);
    func = @(x) myfunc(x,c1, c2);    % function of x alone
    eta3_P(k) = fzero(func,[0.5 1.5]);
    fA3_P(k)  = 1/(1-1/eta3_P(k));
    RPP3(k)   = (1-1/fA3_P(k))*RAP(k);
    DRP3(k)   = RAP(k)-RPP3(k);
    gamma(k)  = DRP3(k)/DRP1(k);
end

Planet =["Jupiter"; "Saturn "; "Uranus "; 'Neptun '];
for k=1:4
fprintf('\n %s  \t %5.2f   \t %5.2f   \t %5.2f  \t %5.2f  \t %5.3f   \t%5.0f km \t%5.0f km \t%5.0f km \t%5.0f km', ...
    Planet(k,:), fA1_P(k),  fA2_P(k), fA3_P(k), fAexp_P(k), gamma(k), DRP1(k)/1000, ...
    DRP2(k)/1000, DRP3(k)/1000, -(RPPex(k)-RAP(k))/1000);
% fprintf("\n");
end

fprintf('\n');
fprintf('\n');
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
