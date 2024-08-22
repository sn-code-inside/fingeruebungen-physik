% -------------------------------------------------------------------------
% FreReturnTrajectoryApollo13.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet für die Parameter der Free Return Trajectorie von Apollo 13.
% 
% -------------------------------------------------------------------------
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
% Ausgangsparameter
h0      = 169609;
r0      = h0+RE;


%AB und Parameter Raumschiff für DGL
P1.RE    = RE;
P1.RM    = RM;
P1.h0    = h0;
P1.muE   = muE;
P1.muM   = muM;
P1.omegaM= omegaM;


%% Berechnung ToF in Patch Conic
% ToF bis zur SOI des Mondes für verschiedene Anfangsparameter 
% Suche nach optimalen Parametern für kleine Periselenium und kleinem
% Perigäum bei Return

%% Berechnung ToF für Patch Conic
% ToF bis zur SOI des Mondes für verschiedene Anfangsparameter 
% Suche nach optimalen Parametern für kleine Periselenium und kleinem
% Perigäum bei Return

%Variiert werden v0 und lambda 1

gamma0  = 0;
lambda1 = deg2rad(linspace(30,90,601));
% Energie Drehimpuls @ TLI

v01      = 10960;
kend = 4;
[TOFh,~,~,hPM,lgdstr]=PatchConic(kend,muE,muM,RSOIM,r0,rEM,...
                                              omegaM,vM,RM,lambda1,v01);

%% Graphik für Ergebnisse aus Patch Conic
figure('name','Patch Conic Näherung')
yyaxis left
hold on
for k=1:kend
    hp(k)=plot(rad2deg(lambda1(:)),hPM(k,:)/1000,...
        'linestyle',Style(k),'linewidth',2);
end
axis([30, 90, -2000 3000])
yticks([-2000 -1000 0 1000 2000 3000])
legend box off;
ylabel('minimale Höhe über Mondoberfläche in km','FontSize',14)
grid on
set(gca,'FontSize',16);

yyaxis right
hold on
for k=1:kend
    plot(rad2deg(lambda1(:)),TOFh(k,:),...
        'linestyle',Style(k),'linewidth',2);
end
hp2=legend(hp,lgdstr,...
    'location','southeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([30, 90, 50 100])
grid on
ylabel('TOF in h','FontSize',14)
xlabel('\lambda_1 in °','FontSize',14)
grid on
set(gca,'FontSize',16);

%% Berechnung Patch Conic für einen Wert von lambda1, v01

lambda1 = deg2rad(49.7);
kend    = 1;
[TOFh,TOFMh,phi0d,hPM,lgdstr]=PatchConic(kend,muE,muM,RSOIM,r0,rEM,...
                                              omegaM,vM,RM,lambda1,v01);

fprintf('\n');
fprintf('Patch Conic Berechnungen\n');
fprintf('\n');
fprintf('h0             :    %6.1f km \n',  h0/1000);
fprintf('v0             :  %8.1f m/s \n',v01);
fprintf('lambda1        :    %7.2f° \n', rad2deg(lambda1));
fprintf('phi0           :    %7.2f° \n', phi0d);
fprintf('ToF zum Mond   :    %6.1f h \n', TOFh+TOFMh/2 );
fprintf('h_min am  Mond :   %+7.1f km \n',hPM/1000);
fprintf('\n');
%% Numerik

% Wir nutzen als Startwerte die Ergebnisse für den Parametersatz 
lambda1 = deg2rad(49.7);
v0      = v01;
alpha0  = pi/4;      %eigentlich irrelevant, da nur Verschiebung
% Wir verkleinern phi0 etwas gegegenüber Patch Conic. 
% Das macht Sinn, da ja die Erdanziehung auch noch
% etwas im Bereich der Mond-SOI wirkt und umgekehrt. Wir setzen als phi0
% als Variable an.
phi0    = deg2rad(129.3253);
beta0   = phi0+alpha0;

%Mond @TLI 
XM0     = rEM*cos(-alpha0);
YM0     = rEM*sin(-alpha0);
%Rakete @TLI
xR0     = r0*cos(-beta0);
yR0     = r0*sin(-beta0);
dXM0    = -omegaM*rEM*sin(-alpha0);
dYM0    = +omegaM*rEM*cos(-alpha0);
dxR0    = -v0*sin(-beta0);
dyR0    = +v0*cos(-beta0);

AB     = [XM0 dXM0 YM0 dYM0 xR0 dxR0 yR0 dyR0];  % AB Mond und Rakete
tv     = linspace(0,192*3600,50001);             % Max Dauer Berechnung 

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
rRE1  = rRE(IminM:length(rRE));
[min_rRE,IminE] = min(rRE1);

fprintf('\n');
fprintf('Numerische Berechnungen\n');
fprintf('\n');
fprintf('h0             :     %6.1f km \n', h0/1000);
fprintf('v0             :   %8.1f m/s \n', v0);
fprintf('phi0           :     %7.2f°  \n', rad2deg(phi0));
fprintf('ToF zum Mond   :     %6.1f h \n', t(IminM)/3600);
fprintf('Gesamtflugzeit :     %6.1f h \n', t(end)/3600);
fprintf('h_min am  Mond :    %+7.1f km \n',(min_rRM-RM)/1000);
fprintf('h_end zur Erde :   %+8.1f km \n', (min_rRE-RE)/1000);
fprintf('\n');

% Mond, Erde Darstellung
[xRE0,yRE0]=Kreis(RE,0,0);
[xRM0,yRM0]=Kreis(RM,XM(1),YM(1));
[xRMend,yRMend]=Kreis(RM,XM(end),YM(end));
[xRMarr,yRMarr]=Kreis(RM,XM(IminM),YM(IminM));
Iacc = 14093;
figure('Name','Numerische Rechnung Free Return Trajectory')
hold on
hp2(4)=plot(xR/1000,yR/1000,LW,2,LC,Colors(4,:));
hp2(2)=plot(XM/1000,YM/1000,LW,2,LS,':',LC,Colors(10,:));
hp2(1)=plot(xRE0/1000,yRE0/1000,LW,1,LC,Colors(3,:));
plot(xRM0/1000,yRM0/1000,LW,1,LC,Colors(10,:));
plot(xRMend/1000,yRMend/1000,LW,1,LC,Colors(10,:));
plot(xRMarr/1000,yRMarr/1000,LW,1,LC,Colors(10,:));
hp2(3) = line([0 XM0]/1000,[0 YM0]/1000);
line([0 XM(end)]/1000,[0 YM(end)]/1000);
line([0 XM(IminM)]/1000,[0 YM(IminM)]/1000);
hp2(5)=plot(xR(Iacc)/1000,yR(Iacc)/1000,'dr',LW,2,'markersize',6);
axis equal
grid on
hp3=legend(hp2(1:5),'Erde @TLI','Mondbahn','Erde @ return', ...
    'FRT Apollo13', 'Ort der Explosion', ...
    'location','southwest'); 
set(hp3,'FontSize',14,'FontWeight','normal'); 
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
        value = (r1 - 6528e3);     % Abbruch bei  h=150 km
        isterminal = 1;            % Stop Integration
        direction  = -1;            % Richtung von oben
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

function [TOFh,TOFMh,phi0d,hPM,lgdstr]=PatchConic(kend,muE,muM,RSOIM,r0,rEM,...
                                             omegaM,vM,RM,lambda1,v01)
  for k=1:kend
    if kend ==1
        v0 = v01;
    else
        v0      = 10950+(k-1)*5;
    end
    lgdstr(k,:) = sprintf('v_0 = %7.3f km/s',v0/1000);
    H0      = v0^2/2-muE/r0;
    L0      = v0*r0;

    % r, v  @ TLI 
    r1      = sqrt(RSOIM^2+rEM^2-2*RSOIM*rEM*cos(lambda1));
    v1      = sqrt(2*(H0+muE./r1));
    
    % FPA @ SOI Mond
    gamma1  = acos(L0./r1./v1);
    gamma1d = rad2deg(gamma1);
    
    % Phasenwinkel Moon @ SOI
    phi1    = asin(RSOIM*sin(lambda1)./r1);
    
    % TOF
    % Ellipsen-Parameter
    
    p0      = L0^2/muE;
    a0      = -muE/H0/2;
    exz0    = sqrt(1-p0/a0);
    cosups0 = 1;
    ups0    = acos(cosups0);
    cosE0   = (exz0+cosups0)/(1+exz0*cosups0);
    E0      = acos(cosE0);
    sinE0   = sin(E0);
    cosups1 = (p0-r1)/exz0./r1;
    ups1    = acos(cosups1);
    cosE1   = (exz0+cosups1)./(1+exz0.*cosups1);
    E1      = acos(cosE1);
    sinE1   = sin(E1);
    
    fac     = sqrt(a0^3/muE);
    TOF     = fac*(E1- E0 -exz0*(sinE1-sinE0));
    TOFh(k,:)    = TOF/3600;
    
    % Phasenwinkel @ TLI
    phi0    = ups1-ups0-phi1-omegaM*TOF;
    phi0d   = rad2deg(phi0);
    
    % Übergang Mond KOS
    r2      = RSOIM;
    v2      = sqrt(v1.^2+vM^2-2*v1.*vM.*cos(gamma1-phi1)); 
    eps2    = asin(vM*cos(lambda1)./v2 - v1.*cos(lambda1+phi1-gamma1)./v2);
    H2      = v2.^2/2-muM./r2;
    % Falls H2 > 0 dann Hyperbelbahn und auch exz2 >1
    L2      = v2.*r2.*sin(eps2);
    p2      = L2.^2/muM;
    exz2    = sqrt(1+2*H2.*L2.^2./muM^2);
    rPM     = p2./(1+exz2);
    vPM     = sqrt(2*(H2+muM./rPM));
    hPM(k,:)= rPM-RM;
   
    EH2      = acosh((1+r2*(exz2-1)./rPM)./exz2);
    a2B      = r2./(exz2.*cosh(EH2)-1);
    cosups2 = a2B.*(exz2-cosh(EH2))./r2;
    ups2    = acos(cosups2);
    % Flugzeit um Mond (innerhalb SOI)
    TOFM  = 2*sqrt(a2B.^3/muM).*(exz2.*sinh(EH2)-EH2);
    TOFMh(k,:) = TOFM/3600;
  end
end

function [xR,yR]=Kreis(R,XM,YM) 
    ups = linspace(0,2*pi,361);
    xR  = R*cos(ups)+XM;
    yR  = R*sin(ups)+YM;
end

% Ende Funktionen
% -------------------------------------------------------------------------
