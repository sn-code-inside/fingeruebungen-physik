% -------------------------------------------------------------------------
% FreReturTrajectory.m
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


%% Berechnung ToF
% ToF bis zur SOI des Mondes für verschiedene Anfangsparameter 
% Suche nach optimalen Parametern für kleine Periselenium und kleinem
% Perigäum bei Return

% Ausgangsparameter
h0      = 169609;
r0      = h0+RE;
%Variiert werden v0 und lambda 1

gamma0  = 0;
lambda1 = deg2rad(linspace(30,90,601));
% Energie Drehimpuls @ TLI

kend = 4;
for k=1:kend
    v0      = 10950+(k-1)*5;
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
    phi1d   = rad2deg(phi1);
    
    % TOF
    % Ellipsen-Parameter
    
    p0      = L0^2/muE;
    a0      = -muE/H0/2;
    exz0    = sqrt(1-p0/a0);
    cosups0 = 1;
    ups0    = acos(cosups0);
    ups0d   = acosd(cosups0);
    cosE0   = (exz0+cosups0)/(1+exz0*cosups0);
    E0      = acos(cosE0);
    E0d     = acosd(cosE0);
    sinE0   = sin(E0);
    
    cosups1 = (p0-r1)/exz0./r1;
    ups1    = acos(cosups1);
    ups1d   = acosd(cosups1);
    cosE1   = (exz0+cosups1)./(1+exz0.*cosups1);
    E1      = acos(cosE1);
    E1d     = acosd(cosE1);
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
    eps2d   = rad2deg(eps2);
    H2      = v2.^2/2-muM./r2;
    % Falls H2 > 0 dann Hyperbelbahn und auch exz2 >1
    L2      = v2.*r2.*sin(eps2);
    p2      = L2.^2/muM;
    exz2    = sqrt(1+2*H2.*L2.^2./muM^2);
    rPM     = p2./(1+exz2);
    vPM     = sqrt(2*(H2+muM./rPM));
    hPM(k,:)     = rPM-RM;
    
    EH2      = acosh((1+r2*(exz2-1)./rPM)./exz2);
    a2B      = r2./(exz2.*cosh(EH2)-1);
    cosups2 = a2B.*(exz2-cosh(EH2))./r2;
    ups2    = acos(cosups2);
    ups2d   = acosd(cosups2);
    
    % Flugzeit um Mond (innerhalb SOI)
    TOFM  = 2*sqrt(a2B.^3/muM).*(exz2.*sinh(EH2)-EH2);
    TOFMh(k,:) = TOFM/3600;
end

% % Graphik für Ergebnisse aus Patch Conic
% figure('name','Patch Conic Näherung')
% yyaxis left
% hold on
% for k=1:kend
%     hp(k)=plot(rad2deg(lambda1(:)),hPM(k,:)/1000,...
%         'linestyle',Style(k),'linewidth',2);
% end
% axis([30, 90, 0 2500])
% yticks([0 500 1000 1500 2000 2500])
% legend box off;
% ylabel('minimale Höhe über Mond in km','FontSize',14)
% grid on
% set(gca,'FontSize',16);
% 
% yyaxis right
% hold on
% for k=1:kend
%     plot(rad2deg(lambda1(:)),TOFh(k,:),...
%         'linestyle',Style(k),'linewidth',2);
% end
% hp2=legend(hp,lgdstr,...
%     'location','southeast'); 
% set(hp2,'FontSize',14,'FontWeight','normal'); 
% axis([30, 90, 50 100])
% grid on
% ylabel('TOF in h','FontSize',14)
% xlabel('\lambda_1 in °','FontSize',14)
% grid on
% set(gca,'FontSize',16);

%% Numerik

% Wir nutzen als Startwerte die Ergebnisse für den Parametersatz 
lambda1 = deg2rad(37.4123);
v0      = 10965.300;
TOF     = 72*3600;
% und verkleinern phi0. Das macht Sinn, da ja die Erdanziehung auch noch
% etwas im Bereich der Mond-SOI wirkt und umgekehrt. Wir setzen als phi0
% als variable an.
alpha0  = omegaM*TOF;
phi0    = deg2rad(130.4269)-deg2rad(2.175);
beta0   = alpha0+phi0;

%AB und Parameter Raumschiff für DGL
P1.RE    = RE;
P1.RM    = RM;
P1.h0    = h0;
P1.muE   = muE;
P1.muM   = muM;
P1.omegaM= omegaM;

%Schwerpunkts-Abstand Erdzentrum
RS0 =  MM*rEM/(ME+MM);

%Koordinaten im Schwerpunktsystem
%Mond @TLI 
XM0     = (rEM-RS0)*cos(-alpha0);
YM0     = (rEM-RS0)*sin(-alpha0);
dXM0    = -omegaM*YM0;
dYM0    = +omegaM*XM0;

%Erde @TLI
XE0     = RS0*cos(pi/2+alpha0);
YE0     = RS0*sin(pi/2+alpha0);
dXE0    = +omegaM*XE0;
dYE0    = -omegaM*YE0;

%Rakete @TLI
xR0     = r0*cos(-beta0)+XE0;
yR0     = r0*sin(-beta0)+YE0;
dxR0    = -v0*sin(-beta0);
dyR0    = +v0*cos(-beta0);

% AB Mond, Rakete, Erde
AB     = [XM0 dXM0 YM0 dYM0 xR0 dxR0 yR0 dyR0 XE0 dXE0 YE0 dYE0];      
iend   = 50001;
tv     = linspace(0,3.0*TOF,iend);             % Max Dauer Berechnung 3*TOF

%% Berechnung 
opts     = odeset('AbsTol',1.e-10,'RelTol',1.e-9);
[t,Y]    = SatellitenPfad(P1, tv, AB);

% Raumschiffbahn und Geschwindigkeit
XM   = Y(:,1);
YM   = Y(:,3);
xR   = Y(:,5);
yR   = Y(:,7);
XE   = Y(:,9);
YE   = Y(:,11);
rRM  = sqrt((xR-XM).^2+(yR-YM).^2);
rRE  = sqrt((xR-XE).^2+(yR-YE).^2);

[min_rRM,IminM] = min(rRM);
rRE1  = rRE(IminM:length(rRE));
[min_rRE,IminE] = min(rRE1);

fprintf('\n');
fprintf('h0             :   %6.1f km \n', h0/1000);
fprintf('v0             :   %8.3f km/s \n', v0/1000);
fprintf('phi0           :   %7.2f ° \n', rad2deg(phi0));
fprintf('lambda1        :   %7.2f ° \n', rad2deg(lambda1));
fprintf('ToF zum Mond   :   %6.1f h \n', t(IminM)/3600);
fprintf('Gesamtflugzeit :   %6.1f h \n', t(end)/3600);
fprintf('h_min am  Mond :  %+7.1f km \n', (min_rRM-RM)/1000);
fprintf('h_end zur Erde : %+8.1f km \n', (min_rRE-RE)/1000);
fprintf('\n');

ups  = linspace(0,2*pi,361);
xRE  = RE*cos(ups)+XE0;
yRE  = RE*sin(ups)+YE0;

xREend  = (RE)*cos(ups)+XE(end);
yREend  = (RE)*sin(ups)+YE(end);

figure('Name','Numerische Rechnung Free Return Trajectory')
hold on
plot(XE/1000,YE/1000,LW,1,LC,Colors(3,:),LS,':');
hp(4)=plot(xR/1000,yR/1000,LW,2);
hp(2)=plot(XM/1000,YM/1000,LW,2,LS,':');
hp(1)=plot(xRE/1000,yRE/1000,LW,1,LC,Colors(3,:));
hp(3)=plot(xREend/1000,yREend/1000,LW,1,LC,Colors(3,:),LS,':');
axis equal
grid on
hp2=legend(hp(1:4),'Erde @TLI','Mondbahn','Erde @ return', ...
    'FRT Apollo13', ...
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
    [tout,Yout,te,ye,ie]=ode45(@(t,Y)MoonTrajectoryCoM(t,Y,P1),tv,AB,opt);
    function [value,isterminal,direction] = events(t,Y)
        rRE1  = sqrt((Y(5)-Y(9)).^2 + (Y(7)-Y(11)).^2);
        value = (rRE1 - 0.98*(P1.RE+P1.h0));     % Abbruch bei h=h0
        isterminal = 1;             % Stop Integration
        direction  = 0;            % Richtung von oben
    end
end

function dY = MoonTrajectoryCoM(t, Y, P1)
    % Y(1,12)=XM,dXM,YM,dYM, x, dx, y, dy, XE,dXE,YE,dYE
    % Es muss ein Spaltenvektor zurÃ¼ckgegeben werden 
    dY  = zeros(12,1); 
    muE = P1.muE;
    muM = P1.muM;
    XM  = Y(1);
    YM  = Y(3);
    xR  = Y(5);
    yR  = Y(7);
    XE  = Y(9);
    YE  = Y(11);
 % Berechnung Mondbahn mit Störung durch Erde
    rRE    =  sqrt((xR-XE).^2+(yR-YE).^2);
    rRM    =  sqrt((xR-XM).^2+(yR-YM).^2);
    rME    =  sqrt((XE-XM).^2+(YE-YM).^2);
    dY(1)  =  Y(2);
    dY(2)  =  -muE*(XM-XE)./rME.^3;    
    dY(3)  =  Y(4);
    dY(4)  =  -muE*(YM-YE)./rME.^3;          
    dY(5)  =  Y(6);
    dY(6)  =  -muE*(xR-XE)./rRE.^3-muM*(xR-XM)/rRM.^3;    
    dY(7)  =  Y(8);
    dY(8)  =  -muE*(yR-YE)./rRE.^3-muM*(yR-YM)/rRM.^3;    
    dY(9)  =  Y(10);
    dY(10)  =  -muM*(XE-XM)./rME.^3;    
    dY(11)  = Y(12);
    dY(12)  = -muM*(YE-YM)./rME.^3;          
end
% Ende Funktionen
% -------------------------------------------------------------------------
