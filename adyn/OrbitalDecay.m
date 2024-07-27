% -------------------------------------------------------------------------
% OrbitalDecay.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Berechnung der Effekte der Reibung auf einen niedrig fliegnden 
% Satelliten (orbital Decay)
%
% Beispiel 1) Tiangong-1   2) ISS (passiv)
% -------------------------------------------------------------------------

%% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
GME = 398600.4415;  % in km^3/s^2
RE  = 6378;         % in km

% Parameter Tiangong-1
rho0      = 6.0e-1;       % Dichte in kg/km^3
H0        = 175;          % Parameter H0 in km
H         = 29.45;         % Parameter H in km
mR        = 8.5e3;        % mR in kg
P1.rho0   = rho0;
P1.RE     = RE;
P1.H0     = H0;
P1.H      = H;
P1.Re     = 122;          % Wiedereintrittssphärenradius


%% Check Dichteabhängigkeit
% h    = linspace(0,500,501);
% rho  = calculate_rho(h, H, H0);
% plot(h,rho);

%% Höhenberechnung

tdays    = 80;
tend     = tdays*24*60*60;

% MATLABs Runge-Kutta ode45 Routine für Apogäum
Aeff      = 63.1e-6;      % Aeff in km^2
kappa     = sqrt(GME)*Aeff/mR;
P1.kappa  = kappa;
h0        = 284;            
AB        = h0;          % AB für DGL Apogäum
opts = odeset('AbsTol',1.e-9,'RelTol',1.e-8,'Events',@Ereignis);
[tA,YA]=ode45(@(t,YA, P1)DGL_OrbitalDecay(t,YA,P1),[0 tend],AB,opts,P1);

% MATLABs Runge-Kutta ode45 Routine für Perigäum
Aeff      = 33.025e-6;      % Aeff in km^2
kappa     = sqrt(GME)*Aeff/mR;
P1.kappa  = kappa;
h0        = 265;            
AB        = h0;          % AB für DGL Perigäum
[tP,YP]=ode45(@(t,YP, P1)DGL_OrbitalDecay(t,YP,P1),[0 tend],AB,opts,P1);

% MATLABs Runge-Kutta ode45 Routine für mittlere Höhe
Aeff      = 46.45e-6;      % Aeff in km^2
kappa     = sqrt(GME)*Aeff/mR;
P1.kappa  = kappa;
h0        = 275;            
AB        = h0;          % AB für DGL Perigäum
[tM,YM]=ode45(@(t,YM, P1)DGL_OrbitalDecay(t,YM,P1),[0 tend],AB,opts,P1);


figure('name','Tiangong-1 Orbital Decay')
tA = tA/60/60/24;
% subplot(2,1,1)
plot(tA,YA,'color',Colors(3,:))
hold on
tP = tP/60/60/24;
plot(tP,YP,'color',Colors(2,:))
tM = tM/60/60/24;
plot(tM,YM,'color',Colors(4,:))
line([0 tdays],[P1.Re P1.Re])
ylabel('Höhe h','FontSize',14)
xlabel('Zeit in Tagen seit 01.02.2018','FontSize',14)
grid on;
h1 = title('Orbital Decay Tiangong-1','FontSize',12);
set(h1,'FontSize',14,'FontWeight','normal'); 
h2=legend('Apogäum','Perigäum','mittlere Höhe','Reentry-Höhe'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
axis([0, tdays, 100 300]);
legend box off;
set(gca,'FontSize',16);


%%

% Parameter ISS 2024

rA1 = 427.5+RE;
rP1 = 412.0+RE;
q1  = rA1/rP1;
e1  = (q1-1)/(q1+1);
a1  = rA1/(1+e1);
h1  = a1-RE;

rA2 = 424.0+RE;
rP2 = 412.5+RE;
q2  = rA2/rP2;
e2  = (q2-1)/(q2+1);
a2  = rA2/(1+e2);
h2  = a2-RE;
Delta_a = a2-a1;        % in km
Delta_t = 11*24*60*60;  % in s

Delta_rt = Delta_a/Delta_t;
H0        = 175;          % Parameter H0 in km
H         = 29.6;         % Parameter H  in km
rhocal = calculate_rho(h1, H, H0);
kappa_eff = -Delta_rt/sqrt(a1)/calculate_rho(h1, H, H0);
mR        = 450e3;        % in kg
Aeff      = kappa_eff*mR/sqrt(GME);     % Aeff in km^2
Aeff_m2   = Aeff*1e6;

P1.rho0   = rho0;
P1.RE     = RE;
P1.H0     = H0;
P1.H      = H;
P1.kappa  = Aeff;
P1.Re     = 122;          % Wiedereintrittssphärenradius
h0        = 427;            
AB        = h0;            % AB für DGL
tdays     = 11;
tend      = tdays*24*60*60;

% MATLABs Runge-Kutta ode45 Routine 
opts = odeset('AbsTol',1.e-9,'RelTol',1.e-8);
[t,Y]=ode45(@(t,Y, P1)DGL_OrbitalDecay(t,Y,P1),[0 tend],AB,opts,P1);

figure('name','ISS')
th = t/60/60/24;
plot(th,Y,'color',Colors(3,:))
ylim([h0-15, h0+1]);
ylabel('Höhe h','FontSize',14)
xlabel('Zeit in Tagen','FontSize',14)
grid on;
hp1 = title('Orbital Decay ISS','FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend('ISS 2014'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
xlim([0, tdays])
ylim([h0-5, h0+1]);
legend box off;
set(gca,'FontSize',16);


%%
% DGL
function dY = DGL_OrbitalDecay(~, Y, P1)
    kap = P1.kappa;
    rho = P1.rho0*exp(-(Y-P1.H0)/P1.H);
    dY = -kap*sqrt(P1.RE+Y)*rho;
end

function [value, isterminal, direction] = Ereignis(t, Y, P1)
    value = Y-P1.Re;
    isterminal = 1;
    direction = 0;
end

function rho = calculate_rho(h, a, b)
    rho  = 0.6*exp(-(h-a)/b); 
end
%% Ende Programm


