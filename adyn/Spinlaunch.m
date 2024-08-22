% -------------------------------------------------------------------------
% Spinlaunch.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm buntersucht die physikalischen Bedingungen für einen
% SpinLaunch
% Alle Rechnungen in kg m s
%
%--------------------------------------------------------------------------


%%
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
LW = 'linewidth';
LS = 'linestyle';
LC = 'color';

%% Initialisierung
% alle Daten in kg, m und s 
G   = 6.67430e-11;                 % G in m^3 /kg /s2
ME  = 5.98e24;                     % Masse der Erde in kg
RE  = 6.378e6;                     % Erdradius in m
gE  = G*ME/RE^2;                   % gE in m/s^2
muE = G*ME;
OmegaE = 2*pi/86165;               % Erdrotation in rad/s

%%-------------------------------------------------------------------------
% Parameterfestlegung
% --------------------------------------------------------------------------

% Parameter der Raketen

% Daten für SpinLaunch
mLA     = 100;                 % Nutzlast in kg
mSA     = 25;                  % Strukturmasse der Stufe in kg
mTA     = 200;                 % Treibstoffsmasse
tBA     = 150;                 % Brenndauer in s
vGA     = 5000;                % Gasausströmgeschwindigkeit 
dmA     = mTA./tBA;            % Massenverlustrate in kg/s
FSA     = dmA*vGA;             % Schub in N
ISpezA  = FSA./dmA/gE;         % Spez. Impuls 
m0      = mSA+mTA+mLA;         % Anfangsmasse
v0      = 8000/3.6;            % Anfangsgeschwindigkeit in m/s
rho0    = 1.752;               % Dichte in kg/m^3
H0      = 6.7e3;               % H0 in km Barometr. Höhenformel
cwD     = 0.10*0.25;           % Widerstandsbeiwert x Fläche in m^2 

% frei wählbare Startparameter zur Bahnoptimierung
tNick    = 60;                 % Zeit bis Einsetzen Schub
Nick_W   = deg2rad(31.5);        % gewählte Ausgangsneigung
fracBurn = 0.50;               % Treibstoffuseanteil Stufe 

Parastr(1,:)= sprintf('Nutzlast : %6.1f kg', mLA);
Parastr(2,:)= sprintf('Nickkorr.: %6.1f ° ', rad2deg(Nick_W));
Parastr(3,:)= sprintf('Used Fuel: %6.3f   ', fracBurn);

% Startpunktgeographische Daten
% Cape Canaveral, 
lat = deg2rad(28.383); 
lon = deg2rad(-80.6);
% Startazimut
azim = deg2rad(90);
R_ax = RE*cos(lat);

% Umfangsgeschwindigkeit auf Erdoberfläche am Startort
vE0x = -OmegaE*R_ax*sin(lon);
vE0y = +OmegaE*R_ax*cos(lon);
vE0z = 0;

% Anfangswerte 
x0          = RE*cos(lat)*cos(lon);
y0          = RE*cos(lat)*sin(lon);
z0          = RE*sin(lat);

% Orientierung der Rakete beim Start (Richtung des Schubs)
vecS(1) = x0/RE;
vecS(2) = y0/RE;
vecS(3) = z0/RE;


% Anfangsgeschwindigkeit (= Mitrotation mit Erde + Anfangswert senkrecht
% nach oben, um numerische Fehler in der DGL zu vermeiden)
vx0 = v0*vecS(1) + vE0x;
vy0 = v0*vecS(2) + vE0y;
vz0 = v0*vecS(3) + vE0z;

% Parameter der DGL
P1.gE       = gE;
P1.muE      = muE;
P1.H0       = H0;
P1.rho0     = rho0;
P1.RE       = RE;
P1.cwD      = cwD;
P1.dmA      = dmA;
P1.tBA      = tBA;
P1.FSA      = FSA;
P1.OmegaE   = OmegaE;
P1.vS(1)    = vE0x;
P1.vS(2)    = vE0y;
P1.vS(3)    = vE0z;

%%--------------------------------------------------------------------------
% Lösen der Differentialgleichung
opts   = odeset('AbsTol',1.e-09,'RelTol',1.e-08);
ispan  = 1001;    % Anzahl der berechneten Zeitschritte pro einzelner Stufe

%--------------------------------------------------------------------------

% Raketenausrichtung:
% Bestimmung der Geschwindigkeit, die sich nicht aus der mitgenommenen
% Rotation der Erde ergibt (= Geschwindigkeitszunahme durch 
% Beschleunigung), denn diese wurde durch den Schub erzeugt und gibt
% die Orientierung der Rakete an.
% vFlugv = [vx1(end)-vE0x; vy1(end)-vE0y; vz1(end)-vE0z];
vFlugv = [vx0; vy0; vz0];
vFlug  = norm(vFlugv);
% Es folgt die Drehung um den Nickkwinkel in Richtung des Startazimuts,
% was mit folgender Drehmatrix RNkorr erreicht wird, wobei erst die
% Nickkorrektur und dann die Ausrichtung auf das Azimut erfolgt
RN = [cos(Nick_W) -sin(Nick_W) 0; sin(Nick_W) cos(Nick_W) 0; 0 0 1];
Ra = [1 0 0; 0 cos(azim) -sin(azim); 0 sin(azim) cos(azim)];
RNkorr = Ra*RN;
% Diese Drehmatrix geht jedoch davon aus, dass die Rakete vor der
% Nickkorrektur entlang der x-Richtung ausgerichtet ist. Sie muss
% in das richtige (sphärische) Koordinatensystem transformiert werden:
% Siehe auch MATLAB File KOS.m
S = [cos(lat)*cos(lon) -sin(lat)*cos(lon) -sin(lon); ...
     cos(lat)*sin(lon) -sin(lat)*sin(lon) +cos(lon); ...
     sin(lat)          +cos(lat)          0];
RNkorrS = S*RNkorr*transpose(S);
% Rakete nach Nickkorrektur
vNick  = RNkorrS*vFlugv;
vNickx = vNick(1) + vE0x;
vNicky = vNick(2) + vE0y;
vNickz = vNick(3) + vE0z;

% Berechnung Anfangsphase 
AB0       = [vNickx;vNicky;vNickz;x0;y0;z0;m0];     % AB für DGL
tspan0    = linspace(0,tNick,ispan);
opts      = odeset('AbsTol',1.e-09,'RelTol',1.e-08,'Events',@gammanull);
[t1,Y,te,ye,ie] = ode45(@(t,Y) DGL_FreeTraj3D(t,Y,P1), tspan0, AB0, opts);
vx1(:)    = Y(:,1);
vy1(:)    = Y(:,2);
vz1(:)    = Y(:,3);
x1(:)     = Y(:,4);
y1(:)     = Y(:,5);
z1(:)     = Y(:,6);
m1(:)     = Y(:,7);
h1(:)     = vecnorm([x1;y1;z1])-RE;
% 


% Berechnung Brennstufe 1 bis fracburn
back      = length(t1)-0;
AB1       = [vx1(end);vy1(end);vz1(end);x1(end);y1(end);z1(end);m1(end)];   
P1.Stage  = 1;
tspan1    = linspace(t1(end),t1(end)+tBA*fracBurn,ispan);
opts      = odeset('AbsTol',1.e-09,'RelTol',1.e-08,'Events', ...
        @(t,Y)absturz(t,Y,P1));
P1.S(:)   = vecS(:);
[t2,Y,te,ye,ie] = ode45(@(t,Y) DGL_Thrust3D(t,Y,P1), tspan1, AB1, opts);
vx2(:)    = Y(:,1);
vy2(:)    = Y(:,2);
vz2(:)    = Y(:,3);
x2(:)     = Y(:,4);
y2(:)     = Y(:,5);
z2(:)     = Y(:,6);
m2(:)     = Y(:,7);
h2(:)     = vecnorm([x2;y2;z2])-RE;

% Berechnung Brennstufe 1 bis Burn-Out und Check auf eventueller Absturz
AB2       = [vx2(end);vy2(end);vz2(end);x2(end);y2(end);z2(end);m2(end)-mSA];   
tspan2    = linspace(t2(end),t1(end)+tBA,ispan);
opts      = odeset('AbsTol',1.e-09,'RelTol',1.e-08,'Events', ...
        @(t,Y)absturz(t,Y,P1));
[t3,Y,te,ye,ie] = ode45(@(t,Y) DGL_Thrust3D(t,Y,P1), tspan2, AB2, opts);
vx3(:)    = Y(:,1);
vy3(:)    = Y(:,2);
vz3(:)    = Y(:,3);
x3(:)     = Y(:,4);
y3(:)     = Y(:,5);
z3(:)     = Y(:,6);
m3(:)     = Y(:,7);
h3(:)     = vecnorm([x3;y3;z3])-RE;

if te >0
    fprintf('\n')
    fprintf('Parameter schlecht gewählt !\n')
    fprintf('Rakete stürzt ab.\n')
    fprintf('\n')
    return
end



% Freie Bewegung
tUmlauf   = 2*pi/sqrt(muE/(RE+200000)^3);
AB3       = [vx3(end);vy3(end);vz3(end);x3(end);y3(end);z3(end);m3(end)];
tspan3    = linspace(t3(end),tUmlauf/2,501);
opts      = odeset('AbsTol',1.e-09,'RelTol',1.e-08,'Events', ...
        @(t,Y)absturz(t,Y,P1));
[t4,Y,te,ye,ie] = ode45(@(t,Y) DGL_FreeTraj3D(t,Y,P1), tspan3, AB3, opts);
vx4(:)    = Y(:,1);
vy4(:)    = Y(:,2);
vz4(:)    = Y(:,3);
x4(:)     = Y(:,4);
y4(:)     = Y(:,5);
z4(:)     = Y(:,6);
m4(:)     = Y(:,7);
h4(:)     = vecnorm([x4;y4;z4])-RE;


%%--------------------------------------------------------------------------
% Graphische Ausgabe
%--------------------------------------------------------------------------


%---------------------------------------------------------------------
% Übersichts-Plots Startphase
taxend = t4(end);  %Zeitdauer für den Plot

figure('name','Übersicht-Plots ohne freie Trajektorie')
subplot(2,2,1)
yyaxis left
hold on
p(1)=plot(t1,h1/1000,LC,Colors(2,:),LW,2,LS,'-');
p(2)=plot(t2,h2/1000,LC,Colors(4,:),LW,2,LS,'-');
p(3)=plot(t3,h3/1000,LC,Colors(3,:),LW,2,LS,'-');
p(4)=plot(t4,h4/1000,LC,Colors(6,:),LW,2,LS,'-');
lgdstr(1,:)='Spinlaunch  Phase    ';
lgdstr(2,:)='Stufe 1 bis FracBurn ';
lgdstr(3,:)='Burn nach Abtrennung ';
lgdstr(4,:)='Free Trajectory      ';
str= "h(t) - , \gamma(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
axis([0, taxend, 0, 250]);
xlabel('t in s','FontSize',14); 
ylabel('h in km','FontSize',14)
grid on
set(gca,'FontSize',16);

yyaxis right
vloc   = sqrt(vx1.^2 + vy1.^2 + vz1.^2);
Rloc   = sqrt(x1.^2 + y1.^2 + z1.^2);
gamma1 = rad2deg(pi/2-acos((x1.*vx1+y1.*vy1+z1.*vz1)./vloc./Rloc));
vloc   = sqrt(vx2.^2 + vy2.^2 + vz2.^2);
Rloc   = sqrt(x2.^2 + y2.^2 + z2.^2);
gamma2 = rad2deg(pi/2-acos((x2.*vx2+y2.*vy2+z2.*vz2)./vloc./Rloc));
vloc   = sqrt(vx3.^2 + vy3.^2 + vz3.^2);
Rloc   = sqrt(x3.^2 + y3.^2 + z3.^2);
gamma3 = rad2deg(pi/2-acos((x3.*vx3+y3.*vy3+z3.*vz3)./vloc./Rloc));
%
vloc   = sqrt(vx4.^2 + vy4.^2 + vz4.^2);
Rloc   = sqrt(x4.^2 + y4.^2 + z4.^2);
gamma4 = rad2deg(pi/2-acos((x4.*vx4+y4.*vy4+z4.*vz4)./vloc./Rloc));
hold on
plot(t1,gamma1,LC,Colors(2,:),LW,2,LS,'--');
plot(t2,gamma2,LC,Colors(4,:),LW,2,LS,'--');
plot(t3,gamma3,LC,Colors(3,:),LW,2,LS,'--');
plot(t4,gamma4,LC,Colors(6,:),LW,2,LS,'--');
line([0,taxend],[0,0]);
axis([0, taxend, -10 90]);
ylabel('\gamma','FontSize',14)
hp2=legend(p(1:4),lgdstr(1:4,:),'location','south'); 
set(hp2,'FontSize',14,'FontWeight','normal');
ylim([-10 90])
legend box off;
grid on
set(gca,'FontSize',16);


subplot(2,2,2)
hold on
vK = sqrt(muE/(h2(end)+RE))/1000;
v1 = vecnorm([vx1;vy1;vz1]);
v2 = vecnorm([vx2;vy2;vz2]);
v3 = vecnorm([vx3;vy3;vz3]);
v4 = vecnorm([vx4;vy4;vz4]);
p(1)=plot(t1,v1/1000,LC,Colors(2,:),LW,2);
p(2)=plot(t2,v2/1000,LC,Colors(4,:),LW,2);
p(3)=plot(t3,v3/1000,LC,Colors(3,:),LW,2);
p(4)=plot(t4,v4/1000,LC,Colors(6,:),LW,2);
p(5)=line([0,taxend],[vK, vK], LW,2, LS,':');
lgdstr(5,:)='Bahngeschw. 200 km   ';
str= "v(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:5),lgdstr(1:5,:),'location','south'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, taxend, 0 , max(10,1.2*v1(end)/1000)]);
xlabel('t in s','FontSize',14); 
ylabel('v in km/s','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);

subplot(2,2,3)
hold on
p(1)=plot(x1/1000,h1/1000,LC,Colors(2,:),LW,2);
p(2)=plot(x2/1000,h2/1000,LC,Colors(4,:),LW,2);
p(3)=plot(x3/1000,h3/1000,LC,Colors(3,:),LW,2);
p(4)=plot(x4/1000,h4/1000,LC,Colors(6,:),LW,2);
hp2=legend(p(1:4),lgdstr(1:4,:),'location','northeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
str= "h(x) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
axis([x1(1)/1000, RE/1000, 0 250]);
xlabel('x in km','FontSize',14); 
ylabel('h in km','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);


subplot(2,2,4)
plot(t1,m1,LC,Colors(2,:),LW,2);
hold on
p(1)=plot(t1,m1,LC,Colors(2,:),LW,2);
p(2)=plot(t2,m2,LC,Colors(4,:),LW,2);
p(3)=plot(t3,m3,LC,Colors(3,:),LW,2);
plot(t4,m4,LC,Colors(6,:),LW,2);
tall = linspace(0,taxend);
p(4)=plot(tall,mLA*tall./tall,LC,Colors(8,:),LW,2);
str= "m(t) - Graphik";
hp1 = title(str,'FontSize',12);
lgdstr(4,:)='Nutzlast             ';
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:4),lgdstr(1:4,:),'location','northeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, tNick+2*tBA, 0.8*mLA 1.2*m0]);
xlabel('t in s','FontSize',14); 
ylabel('m in kg','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);

%%--------------------------------------------------------------------------
%  Funktionen
%--------------------------------------------------------------------------

function dY = DGL_Thrust3Dfixed(t, Y, P1)
% Diese Variante geht von einer strikt senkrecht nach oben ausgerichteten
% Rakete aus. Für den Start ist es notwendig, erst Geschwindigkeit
% senkrecht nach oben zu gewinnen, um nicht abzudriften und um die 
% aerodynamischen Belastungen auf die Rakete zu minimieren.
% Y(1):vx, Y(2):vy, Y(3):vz, Y(4):x Y(5):y Y(6):z Y(7):m
R      = norm([Y(4) Y(5) Y(6)]);
rho    = P1.rho0*exp(-(R-P1.RE)/P1.H0);
% Geschwindigkeit relativ zur Atmosphäre für Reibung
vE_rot = P1.OmegaE*sqrt(Y(4)^2+Y(5)^2);
vRv    = [Y(1)+vE_rot*sin(atan2(Y(5),Y(4))) ...
         Y(2)-vE_rot*cos(atan2(Y(5),Y(4))) Y(3)];
vR     = norm(vRv);
F_D    = 0.5*rho*P1.cwD*vR^2;
C1     = P1.FSA/Y(7);
C2     = -F_D/vR/Y(7);
C3     = -P1.muE/R^3;
dY     = [C1*P1.S(1)+C2*vRv(1)+C3*Y(4);        
          C1*P1.S(2)+C2*vRv(2)+C3*Y(5);  
          C1*P1.S(3)+C2*vRv(3)+C3*Y(6);  
          Y(1);
          Y(2);
          Y(3);
          -P1.dmA];
end

function dY = DGL_Thrust3D(t, Y, P1)
% Diese Variante der Differentialgleichung geht davon aus, dass
% die Rakete sich entlang ihrer Geschwindigkeit ausrichtet, wobei
% die Geschwindigkeit, die sie aufgrund der Erdrotation mitnimmt
% nicht berücksichtigt wird.
% Y(1):vx, Y(2):vy, Y(3):vz, Y(4):x Y(5):y Y(6):z Y(7):m
R      = norm([Y(4) Y(5) Y(6)]);
rho    = P1.rho0*exp(-(R-P1.RE)/P1.H0);
% Geschwindigkeit relativ zur Atmosphäre für Reibung
vE_rot = P1.OmegaE*sqrt(Y(4)^2+Y(5)^2);
vRv    = [Y(1)+vE_rot*sin(atan2(Y(5),Y(4))) ...
          Y(2)-vE_rot*cos(atan2(Y(5),Y(4))) Y(3)];
vR     = norm(vRv);
% Geschwindigkeit, abzüglich der initialen Geschwindigkeit aufgrund
% der Erdrotation
vRSv   = [Y(1)-P1.vS(1) Y(2)-P1.vS(2) Y(3)-P1.vS(3)];
vRS    = norm(vRSv);
F_D    = 0.5*rho*P1.cwD*vR^2;
C1     = P1.FSA/vRS/Y(7);
C2     = -F_D/vR/Y(7);
C3     = -P1.muE/R^3;
dY     = [C1*vRSv(1)+C2*vRv(1)+C3*Y(4);        
          C1*vRSv(2)+C2*vRv(2)+C3*Y(5);  
          C1*vRSv(3)+C2*vRv(3)+C3*Y(6);  
          Y(1);
          Y(2);
          Y(3);
          -P1.dmA];
end

function dY = DGL_FreeTraj3D(t, Y, P1)
% Y(1):vx, Y(2):vy, Y(3):vz, Y(4):x Y(5):y Y(6):z Y(7):m
R      = norm([Y(4) Y(5) Y(6)]);
rho    = P1.rho0*exp(-(R-P1.RE)/P1.H0);
% Geschwindigkeit relativ zur Atmosphäre für Reibung
vE_rot = P1.OmegaE*sqrt(Y(4)^2+Y(5)^2);
vRv    = [Y(1)+vE_rot*sin(atan2(Y(5),Y(4))) ...
          Y(2)-vE_rot*cos(atan2(Y(5),Y(4))) Y(3)];
vR     = norm(vRv);
F_D    = 0.5*rho*P1.cwD*vR^2;
C1     = -F_D/vR/Y(7);
C2     = P1.muE/R^3;
dY     = [C1*vRv(1)-C2*Y(4);        
          C1*vRv(2)-C2*Y(5);  
          C1*vRv(3)-C2*Y(6);  
          Y(1);
          Y(2);
          Y(3);
          0];
end

% Halte an, wenn gamma = 0 erreicht wird
function [value,isterminal,direction] = gammanull(t, Y)
        vloc  = sqrt(Y(1)^2 + Y(2)^2 + Y(3)^2);
        Rloc  = sqrt(Y(4)^2 + Y(5)^2 + Y(6)^2);
        value = pi/2-acos((Y(1)*Y(4)+Y(2)*Y(5)+Y(3)*Y(6)) ...
            ./vloc./Rloc);  % identifiziere gamma = 0
        isterminal = 1;     % Integration anhalten
        direction = 0;      % gamma = 0 wird von oben erreicht
end

% Halte an, wenn h = 0 erreicht wird
function [value,isterminal,direction] = absturz(t, Y, P1)
        value = sqrt(Y(4)^2 + Y(5)^2 + Y(6)^2) - P1.RE;  % identifiziere h=0
        isterminal = 1;     % Integration anhalten
        direction = 0;      % h = 0 wird von oben erreicht
end

