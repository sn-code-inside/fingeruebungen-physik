% -------------------------------------------------------------------------
% FreReturTrajectory.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet für verschiedene Parameter die Free Return Trajectorie um den
% Mond.
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
% alle Daten in kg, m und s 
G   = 6.67430e-11;                 % G in m^3 /kg /s^2
ME  = 5.977e24;                    % Masse der Erde in kg
MM  = 7.348e22;                    % Masse Mond in kg
RE  = 6.378e6;                     % Erdradius in m
RM  = 1.7374e6;                    % Mondradius in m
rEM = 384400e3;                    % Abstand Erde-Mond in m
omegaM = 2.649e-6;                 % Umlauffrequenz Mond um Erde in rad/s
vM     = omegaM*rEM;               % Bahngeschwindigkeit Mond in m/s
RSOIM  = rEM*(MM/ME)^(2/5);        % SOI-Radius Mond
muE    = G*ME;
muM    = G*MM;


%% Berechnung ToF
% ToF bis zur SOI des Mondes für verschiedene Anfangsparameter 

% Ausgangsparameter
h0      = 317000;
r0      = h0+RE;
v0      = 10.83697*1000;
gamma0  = 0;
lambda1 = deg2rad(45.0);
% Energie Drehimpuls @ TLI
H0      = v0^2/2-muE/r0;
L0      = v0*r0;
% r, v  @ TLI 
r1      = sqrt(RSOIM^2+rEM^2-2*RSOIM*rEM*cos(lambda1));
v1      = sqrt(2*(H0+muE/r1));
% FPA @ SOI Mond
gamma1  = acos(L0/r1/v1); gamma1d = rad2deg(gamma1);

% Phasenwinkel Moon @ SOI
phi1    = asin(RSOIM*sin(lambda1)/r1); phi1d   = rad2deg(phi1);

% TOF
% Ellipsen-Parameter
p0      = L0.^2/muE;
a0      = -muE./H0/2;
exz0    = sqrt(1-p0./a0);
cosups0 = 1;
ups0    = acos(cosups0);
ups0d   = acosd(cosups0);
cosE0   = (exz0+cosups0)./(1+exz0.*cosups0);
E0      = acos(cosE0);
E0d     = acosd(cosE0);
sinE0   = sin(E0);

cosups1 = (p0-r1)./exz0/r1;
ups1    = acos(cosups1);
ups1d   = acosd(cosups1);
cosE1   = (exz0+cosups1)./(1+exz0.*cosups1);
E1      = acos(cosE1);
E1d     = acosd(cosE1);
sinE1   = sin(E1);

fac     = sqrt(a0.^3/muE);
TOF     = fac*(E1- E0 -exz0*(sinE1-sinE0));
TOFh    = TOF/3600;

% Phasenwinkel @ TLI
phi0    = ups1-ups0-phi1-omegaM*TOF; phi0d   = rad2deg(phi0);



%AB und Parameter Raumschiff
alpha0 = omegaM*TOF;
phi0   = phi0*0.98646;
beta0  = alpha0+phi0;
P1.RE    = RE;
P1.h0    = h0;
P1.muE   = muE;
P1.muM   = muM;
P1.omegaM= omegaM;


X0     = rEM*cos(-alpha0);
Y0     = rEM*sin(-alpha0);
x0     = r0*cos(-beta0);
y0     = r0*sin(-beta0);
dX0    = -omegaM*rEM*sin(-alpha0);
dY0    = +omegaM*rEM*cos(-alpha0);
dx0    = -v0*sin(-beta0);
dy0    = +v0*cos(-beta0);

AB     = [X0 dX0 Y0 dY0 x0 dx0 y0 dy0];       % AB Mond und Rakete
tv     = linspace(0,3.5*TOF,10001);             % Geschätze Dauer Vorbeiflug

%% Berechnung 
opts     = odeset('AbsTol',1.e-9,'RelTol',1.e-8);
[t,Y]    = SatellitenPfad(P1, tv, AB);

% Raumschiffbahn und Geschwindigkeit
XM   = Y(:,1);
YM   = Y(:,3);
xR   = Y(:,5);
yR   = Y(:,7);
rRM  = sqrt((xR-XM).^2+(yR-YM).^2);
rRE  = sqrt((xR).^2+(yR).^2);

[min_rRM,IminM] = min(rRM);
[min_rRE,IminE] = min(rRE);

fprintf('\n');
fprintf('h0             :  %6.1f km \n', h0/1000);
fprintf('v0             :  %8.3f km/s \n', v0/1000);
fprintf('phi0           :  %7.2f ° \n', rad2deg(phi0));
fprintf('lambda1        :  %7.2f ° \n', rad2deg(lambda1));
fprintf('ToF zum Mond   :  %6.1f h \n', t(IminM)/3600);
fprintf('Gesamtflugzeit :  %6.1f h \n', t(end)/3600);
fprintf('h_min am  Mond : %+7.1f km \n', (min_rRM-RM)/1000);
fprintf('h_min zur Erde : %+7.1f km \n', (min_rRE-RE)/1000);
fprintf('\n');

ups  = linspace(0,2*pi,361);
xRE  = RE*cos(ups);
yRE  = RE*sin(ups);

figure('Name','Free Return Trajectory')
hold on
hp(1)=plot(xR/1000,yR/1000);
hp(2)=plot(xRE/1000,yRE/1000);
hp(3)=plot(XM/1000,YM/1000);
hp(4)=plot(XM/1000,YM/1000);
axis equal
grid on
axis equal
hp2=legend(hp,'Free Return Trajectory','Erde','Mondbahn','SOI Mond',...
    'location','southwest'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
legend box off;
xlabel('x in km','FontSize',14); 
ylabel('y in km','FontSize',14)
set(gca,'FontSize',16);


% -------------------------------------------------------------------------
% Funktionen
% DGL
function [tout,Yout] = SatellitenPfad(P1, tv, AB)
    opt = odeset('AbsTol',1.e-9,'RelTol',1.e-8,'Events',@events);
    [tout,Yout,te,ye,ie]=ode45(@(t,Y)MoonTrajectory(t,Y,P1),tv,AB,opt);
    function [value,isterminal,direction] = events(t,Y)
        r1  = sqrt((Y(5)).^2 + (Y(7)).^2);
        value = (r1 - 0.98*(P1.RE+P1.h0));     % Abbruch bei h=h0
        isterminal = 1;             % Stop Integration
        direction  = 1;             % Egal welche Richtung
    end
end

function dY = MoonTrajectory(t, Y, P1)
    % Y(1,8)=XM,dXM,YM,dYM, x, dx, y, dy
    % Es muss ein Spaltenvektor zurÃ¼ckgegeben werden 
    dY     =  zeros(8,1); 
    omegaM = P1.omegaM;
    muE    = P1.muE;
    muM    = P1.muM;
    XM  = Y(1);
    YM  = Y(3);
    xR  = Y(5);
    yR  = Y(7);
    %Vereinfachte Berechnung Mondbahn (=Kreisbahn)
    dY(1)  =  Y(2);
    dY(2)  = -omegaM^2*XM;    
    dY(3)  =  Y(4);
    dY(4)  = -omegaM^2*YM;       
    rR     =  sqrt(xR.^2 + yR.^2);
    rM     =  sqrt((xR-XM).^2+(yR-YM).^2);
    dY(5)  =  Y(6);
    dY(6)  = -muE*xR./rR.^3+muM*(XM-xR)/rM.^3;    
    dY(7)  =  Y(8);
    dY(8)  = -muE*yR./rR.^3+muM*(YM-yR)/rM.^3;  
end
% Ende Funktionen
% -------------------------------------------------------------------------
