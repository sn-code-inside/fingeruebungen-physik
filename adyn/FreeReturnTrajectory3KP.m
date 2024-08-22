%%-------------------------------------------------------------------------
% FreReturnTrajectory3KP.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Himmelsmechanik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Berechnet für die Parameter der Free Return Trajectorie von Apollo 13 
% 
% -------------------------------------------------------------------------

%% Initialisierung

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", "-",":"];
LW = 'linewidth';
LC = 'color';
LS = 'linestyle';

%% Initialisierung
% alle Werte in kg, m und s 
G   = 6.67430e-11;                 % G in m^3 /kg /s^2
ME  = 5.977e24;                    % Masse der Erde in kg
MM  = 7.348e22;                    % Masse Mond in kg
RE  = 6.378e6;                     % Erdradius in m
RM  = 1.7374e6;                    % Mondradius in m
rEM = 384400e3;                    % Abstand Erde-Mond in m
omegaM = 2.649e-6;                 % Umlauffrequenz Mond um Erde in rad/s
muE    = G*ME;
muM    = G*MM;
% Ausgangsparameter
h0      = 169609;
r0      = h0+RE;
beta0   = deg2rad(47.45);
phi0    = pi-beta0;
v0      = 10953.8055;

meffE   = ME/(ME+MM);
meffM   = MM/(ME+MM);

%Schwerpunkts-Abstand Erdzentrum
XE0 =  -meffM*rEM;
XM0 =  rEM*(1-meffM);

P1.muE  = muE;
P1.muM  = muM;
P1.rEM  = rEM;
P1.meffE  = meffE;
P1.meffM  = meffM;
P1.omegaM = omegaM;
P1.XE0  = XE0;
P1.XM0  = XM0;
P1.RE   = RE;


%% Anfangspositionen und -geschwindigkeitender Körper

xR0     = r0*cos(-phi0)+XE0;
yR0     = r0*sin(-phi0);
vxR0    = -v0*sin(-phi0);
vyR0    = v0*cos(-phi0);

AB   = [xR0, vxR0, yR0, vyR0];
tend = 192*3600;
tv   = linspace(0,tend,10001);

%% Berechnung Trajektorien über ode45
opts = odeset('AbsTol',1.e-09,'RelTol',1.e-08);
[tA,YA]    = SatellitenPfad(P1, tv, AB);

% [tA,YA]=ode45(@(t,YA, P1)SatellitenPfad(P1, tv, AB),[0 tend],AB,opts,P1);

%% Graphische Ausgabe

xR   = YA(:,1);
yR   = YA(:,3);

rRM  = sqrt((xR-XM0).^2+(yR-0).^2);
rRE  = sqrt((xR-XE0).^2+(yR-0).^2);

[min_rRM,IminM] = min(rRM);
rRE1  = rRE(IminM:length(rRE));
[min_rRE,IminE] = min(rRE1);

fprintf('\n');
fprintf('Numerische Berechnungen\n');
fprintf('\n');
fprintf('h0             :     %6.1f km \n', h0/1000);
fprintf('v0             :   %8.1f m/s \n', v0);
fprintf('phi0           :     %7.2f° \n', rad2deg(phi0));
% fprintf('ToF zum Mond   :     %6.1f h \n', t(IminM)/3600);
% fprintf('Gesamtflugzeit :     %6.1f h \n', t(end)/3600);
fprintf('h_min am  Mond :    %+7.1f km \n',(min_rRM-RM)/1000);
fprintf('h_end zur Erde :   %+8.1f km \n', (min_rRE-RE)/1000);
fprintf('\n');


[xRE0,yRE0]=Kreis(RE,XE0,0);
[xRM0,yRM0]=Kreis(RM,XM0,0);


figure('name','Trajektorie')
hold on
hp(1)=plot(xRE0/1000,yRE0/1000,LW,2,LC,Colors(3,:));
hp(2)=plot(xRM0/1000,yRM0/1000,LW,2,LC,Colors(10,:));
hp(3)=plot(xR/1000,yR/1000,LW,2,LC,Colors(4,:));
axis([-10*RE rEM+10*RM 1.05*min(yR) 1.05*max(yR)]/1000)
axis equal
grid on
xlabel('x in km')
ylabel('y in km')
set(gca,'FontSize',14);

f = text(0, 10*RE/1000, 't = 0 s');
h = plot(xR(1)/1000,yR(1)/1000,'d','color', Colors(4,:));
explosion = false;
for k=1:10:length(xR)
  if tA(k)/3600 > 56 && ~explosion 
     explosion = true;
     hp(4)=plot(xR(k)/1000,yR(k)/1000,'dk',LW,2);
  end
  h.Visible = 'off';
  f.Visible = 'off';
  strZeit = sprintf('t = %10.2f h',tA(k)/3600);
  h = plot(xR(k)/1000,yR(k)/1000,'d',LC, Colors(4,:));%Position Rakete
  f = text(0, 10*RE/1000, strZeit);
  h.Visible = 'on';
  f.Visible = 'on';
  pause(0.05)
end
legend(hp, 'Erde', 'Mond', 'Free Return Trajectory','Ort der Explosion',...
    'location', 'northeast')
legend box off
grid on
set(gca,'FontSize',14);


%% DGL -------------------------------------------------------------------

% DGL
function [tout,Yout] = SatellitenPfad(P1, tv, AB)
    opt = odeset('AbsTol',1.e-9,'RelTol',1.e-8,'Events',@events);
    [tout,Yout,~,~,~]=ode45(@(t,Y)DGL_3BodySystem(t,Y,P1),tv,AB,opt);
    function [value,isterminal,direction] = events(~,Y)
        rRE1  = sqrt((Y(1)-P1.XE0).^2 + Y(3).^2);
        value = (rRE1 - 122000-P1.RE);    % Abbruch bei h=122 km
        isterminal = 1;             % Stop Integration
        direction  = 0;             % Richtung von oben
    end
end

function dY = DGL_3BodySystem(~, Y, P1)
    muE = P1.muE;
    muM = P1.muM;
    rEM = P1.rEM;
    meffE = P1.meffE;
    meffM = P1.meffM;
    omegaM= P1.omegaM;
    xR    = Y(1);
    yR    = Y(3);
    vxR   = Y(2);
    vyR   = Y(4);
    rRE  = norm([xR,yR]-[P1.XE0,0]);
    rRM  = norm([xR,yR]-[P1.XM0,0]);
    dY =  [vxR;...
          +2*omegaM*vyR+omegaM^2*xR-muE*(xR+meffM*rEM)/rRE^3-...
                                   muM*(xR-meffE*rEM)/rRM^3;
           vyR;                       
          -2*omegaM*vxR+omegaM^2*yR-muE*yR/rRE^3-...
                                   muM*yR/rRM^3;
         ];
end

function [xR,yR]=Kreis(R,XM,YM) 
    ups = linspace(0,2*pi,361);
    xR  = R*cos(ups)+XM;
    yR  = R*sin(ups)+YM;
end


%% Ende Programm


