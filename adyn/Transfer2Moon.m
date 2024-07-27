% -------------------------------------------------------------------------
% Transfer2Moon.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet verschiedene Parameter u. Trajektorien für einen Flug zum Mond
% 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Initialisierung
% alle Daten in kg, m und s 
G   = 6.67430e-11;                 % G in m^3 /kg /s2
ME  = 5.977e24;          
% Masse der Erde in kg
RE  = 6.378e6;                     % Erdradius in m
muE = G*ME;
rEM = 384400e3;                    % Abstand Erde-Modn in m
omegaM = 2.649e-6;                 % Umlauffrequenz Mond um Erde in rad/s
MM     = 7.3483e22 ;               % Masse Mond in kg
RSOIM  = rEM*(MM/ME)^(2/5);        % SOI-Radius Mond


%% Berechnung ToF
% ToF bis zur SOI des Mondes für verschiedene Anfangsparameter 

h0 = 300;  % Anfangshöhe in km
r0 = RE+h0*1000;

v0 = linspace(10.8,11.2,1001)*1000; % Anfangsgeschwindigkeit in m/s
gamma0  = 0;
lambda1 = deg2rad(linspace(0,90,4));

H0 = v0.^2/2-muE./r0;
L0 = r0*v0*cos(gamma0);



% Ellipsen-Parameter

p    = L0.^2/muE;
a    = -muE./H0/2;
exz  = sqrt(1-p./a);
% e0 = sqrt(1+L0.^2*2.*H0/muE^2);

r1  = sqrt(rEM^2 + RSOIM^2-2*rEM*RSOIM.*cos(lambda1));

cosups0 = (p-r0)./exz/r0;
cosE0 = (exz+cosups0)./(1+exz.*cosups0);
E0    = acos(cosE0);
sinE0 = sin(E0);
for k=1:length(lambda1)
    cosups1(k,:) = (p(:)-r1(k))./exz(:)/r1(k);
end
cosE1 = zeros(length(lambda1), length(v0));
for k=1:length(lambda1)
    for m=1:length(v0)
        cosE1(k,m)= (exz(m)+cosups1(k,m))./(1+exz(m).*cosups1(k,m));
        E1(k,m)   = acos(cosE1(k,m));
        sinE1(k,m)= sin(E1(k,m));
    end
end

% ToF-berechnung
fac = sqrt(a.^3/muE);
for k=1:length(lambda1)
    for m=1:length(v0)
        TOF(k,m)= fac(m).*(E1(k,m)- E0(m)-exz(m)*(sinE1(k,m)-sinE0(m)));
        if ~ isreal(TOF(k,m))
            TOF(k,m) =  NaN;
        end
    end
end

figure('name','TOF');
hold on
for k=1:length(lambda1)
  hp(k)=plot(v0(:)/1000, TOF(k,:)/3600,'linewidth',2);
  lgdstr(k,:) = sprintf(' %4.1f°',rad2deg(lambda1(k)));
end
grid on
str= "Time of Flight zum Mond";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(hp,lgdstr,'location','northeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
legend box off;
axis([min(v0)/1000, max(v0)/1000, 20 100]);
xlabel('v in km/s','FontSize',14); 
ylabel('t in h','FontSize',14)
grid on
set(gca,'FontSize',16);

figure('name','vmin  und TOFmax');
h0     = linspace(0,10000,101)*1000;
r0     = RE+h0;
a0     = (r0+rEM)/2;
vmin   = sqrt(2*muE*(1-1./(1+rEM./r0))./r0);
TOFmin = pi*sqrt(a0.^3/muE);
yyaxis left
hold on
plot((r0(:)-RE)/1000, vmin/1000,'linewidth',2);
axis([0, max(h0)/1000, round(min(vmin)/1000-1) round(max(vmin)/1000+1)]);
grid on
xlabel('h_0 in km','FontSize',14); 
ylabel('v_{min} in km/s','FontSize',14)
yyaxis right
hold on
plot((r0(:)-RE)/1000, TOFmin/3600,'linewidth',2);
ylabel('ToF_{max} in h','FontSize',14)
str= "Time of Flight zum Mond";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);


