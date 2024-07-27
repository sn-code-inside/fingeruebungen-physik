% -------------------------------------------------------------------------
% RakentenStart2Orbit04.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet Start einer dreistufigen Rakete von der Erdoberfläche
% unter Berücksichtigung von Reibung und Nickwinkel-Korrekturen im
% Inertialsystem, das seinen Ursprung im Erdmittelpunkt hat und
% von den Einheitsvektoren der Kugelkoordinaten (Orientiert am
% Startort zum Startzeitpunkt) aufgespannt wird.
% Variabel: Raketenparameter, Nutzlast, Resttreibstoff
% Bsp. Ariane 1
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

% Daten Ariane 1 
mLA     = 1500;
mSA     = [153.870 36.271 9.369]*1e3;
mTA     = [140	33.03 8.238]*1e3;
m0      = sum(mSA)+mLA;
m0A     = [m0	m0-mSA(1)	m0-mSA(1)-mSA(2)];
tBA     = [138	130	563];
FSA     = [2.560 0.71 0.059]*1e6;
dmA     = mTA./tBA;
ISpezA  = FSA./dmA/gE;
mLA     = m0A(3)-mSA(3);
sigmaL  = mLA/m0;

% frei wählbare Startparameter zur Bahnoptimierung
tNick    = 20;                                % Zeit der Nickkorrektur in s
Nick_W_s = deg2rad([1.0,1.5,1.8,1.85,1.9,2.1]); % Mögliche Nickkorrekturen
Nick_W   = Nick_W_s(4);                       % gewählte Nickkorrektur
fracBurn = 1;                                 % Treibstoffuseanteil Stufe 3

Parastr(1,:)= sprintf('Nutzlast : %6.1f kg', mLA);
Parastr(2,:)= sprintf('Nickkorr.: %6.1f ° ', rad2deg(Nick_W));
Parastr(3,:)= sprintf('Used Fuel: %6.3f   ', fracBurn);

% Startpunktgeographische Daten
% 1 - Cape Canaveral, 2 - Baikonur, 3 - Kourou
lats = deg2rad([28.383;45.617;5.158;0]);
lons = deg2rad([-80.6;63.317;-52.65;0]);
% Startazimut
% 1 - Cape Canaveral, 2 - Baikonur, 3 - Kourou
azims = deg2rad([120;80;60]);

% Anfangswerte (Cape Canaveral)
lat  = lats(1); 
lon  = lons(1);
azim = azims(1);
R_ax = RE*cos(lat);

% Umfangsgeschwindigkeit auf Erdoberfläche am Startort
vE0x = -OmegaE*R_ax*sin(lon);
vE0y = +OmegaE*R_ax*cos(lon);
vE0z = 0;

% Anfangswerte 
x0          = RE*cos(lat)*cos(lon);
y0          = RE*cos(lat)*sin(lon);
z0          = RE*sin(lat);
m0          = m0A(1);

% Orientierung der Rakete beim Start (Richtung des Schubs)
vecS(1) = x0/RE;
vecS(2) = y0/RE;
vecS(3) = z0/RE;

% Anfangsgeschwindigkeit (= Mitrotation mit Erde + Anfangswert senkrecht
% nach oben, um numerische Fehler in der DGL zu vermeiden)
vx0 = 0.001*vecS(1) + vE0x;
vy0 = 0.001*vecS(2) + vE0y;
vz0 = 0.001*vecS(3) + vE0z;

% Parameter der DGL
rho0       = 1.752;               % Dichte in kg/m^3
H0         = 6.7e3;               % H0 in km Barometr. Höhenformel
cwD        = 1.5;                 % Widerstandsbeiwert x Fläche in m^2 
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
%--------------------------------------------------------------------------

% MATLABs Runge-Kutta ode45 Routine 
opts   = odeset('AbsTol',1.e-09,'RelTol',1.e-08);
ispan  = 501;    % Anuzal der berechneten Zeitschritte pro einzele Stufe

% Berechnung Anfangsphase (bis Nickkorrektur)
AB0       = [vx0;vy0;vz0;x0;y0;z0;m0];     % AB für DGL
P1.Stage  = 1;
P1.S(:)   = vecS(:);
tspan1    = linspace(0,tNick,101);
[t1,Y]    = ode45(@DGL_Thrust3Dfixed,tspan1, AB0,opts,P1);
vx1(:)    = Y(:,1);
vy1(:)    = Y(:,2);
vz1(:)    = Y(:,3);
x1(:)     = Y(:,4);
y1(:)     = Y(:,5);
z1(:)     = Y(:,6);
m1(:)     = Y(:,7);
h1(:)     = vecnorm([x1;y1;z1])-RE;

% Nickkorrektur:
% Bestimmung der Geschwindigkeit, die sich nicht aus der mitgenommenen
% Rotation der Erde ergibt (= Geschwindigkeitszunahme durch 
% Beschleunigung), denn diese wurde durch den Schub erzeugt und gibt
% die Orientierung der Rakete an.
vFlugv = [vx1(end)-vE0x; vy1(end)-vE0y; vz1(end)-vE0z];
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

% Berechnung erste Stufe bis Burn-Out Stufe 1
% AB für DGL
AB1       = [vNickx;vNicky;vNickz;x1(end);y1(end);z1(end);m1(end)];   
P1.Stage  = 1;
tspan1    = linspace(tNick,tBA(1),ispan);
[t2,Y]    = ode45(@DGL_Thrust3D,tspan1, AB1,opts,P1);
vx2(:)    = Y(:,1);
vy2(:)    = Y(:,2);
vz2(:)    = Y(:,3);
x2(:)     = Y(:,4);
y2(:)     = Y(:,5);
z2(:)     = Y(:,6);
m2(:)     = Y(:,7);
h2(:)     = vecnorm([x2;y2;z2])-RE;

% Berechnung zweite Stufe bis Burn-Out Stufe 2 und Check evtl. Absturz
AB2       = [vx2(end);vy2(end);vz2(end);x2(end);y2(end);z2(end);m0A(2)];   
P1.Stage  = 2;
tspan2    = linspace(t2(end),t2(end)+tBA(2),ispan);
opts      = odeset('AbsTol',1.e-09,'RelTol',1.e-08,'Events', ...
        @(t3,Y)absturz(t3,Y,P1));
[t3,Y,te,ye,ie] = ode45(@(t3,Y) DGL_Thrust3D(t3,Y,P1), tspan2, AB2, opts);
vx3(:)    = Y(:,1);
vy3(:)    = Y(:,2);
vz3(:)    = Y(:,3);
x3(:)     = Y(:,4);
y3(:)     = Y(:,5);
z3(:)     = Y(:,6);
m3(:)     = Y(:,7);
h3(:)     = vecnorm([x3;y3;z3])-RE;

% Abfrage ob bis zum Burn-Out Stufe 2 bzw 3 Absturz erfolgt
% werden kann

if te >0
    fprintf('\n')
    fprintf('Parameter schlecht gewählt !\n')
    fprintf('Rakete stürzt ab.\n')
    fprintf('\n')
    return
end

% Berechnung dritte Stufe bis Burn-Out und Check auf eventueller Absturz
AB3       = [vx3(end);vy3(end);vz3(end);x3(end);y3(end);z3(end);m0A(3)];   
P1.Stage  = 3;
tspan3    = linspace(t3(end),t3(end)+tBA(3)*fracBurn,ispan);
opts      = odeset('AbsTol',1.e-09,'RelTol',1.e-08,'Events', ...
        @(t3,Y)absturz(t3,Y,P1));
[t4,Y,te,ye,ie] = ode45(@(t,Y) DGL_Thrust3D(t,Y,P1), tspan3, AB3, opts);

if te >0
    fprintf('\n')
    fprintf('Parameter schlecht gewählt !\n')
    fprintf('Rakete stürzt ab.\n')
    fprintf('\n')
    return
end

% Berechnung dritte Stufe bis Burn-Out Stufe 3 und FPA 0
AB3       = [vx3(end);vy3(end);vz3(end);x3(end);y3(end);z3(end);m0A(3)];   
P1.Stage  = 3;
tspan3    = linspace(t3(end),t3(end)+tBA(3)*fracBurn,ispan);
opts      = odeset('AbsTol',1.e-09,'RelTol',1.e-08,'Events',@gammanull);
[t4,Y,te,ye,ie] = ode45(@(t,Y) DGL_Thrust3D(t,Y,P1), tspan3, AB3, opts);
vx4(:)    = Y(:,1);
vy4(:)    = Y(:,2);
vz4(:)    = Y(:,3);
x4(:)     = Y(:,4);
y4(:)     = Y(:,5);
z4(:)     = Y(:,6);
m4(:)     = Y(:,7);
h4(:)     = vecnorm([x4;y4;z4])-RE;


% Abfrage ob bis zum Burn-Out Stufe 3 FPA=0 und damit Kreisbahn erreicht
% werden kann
if length(t4) == ispan
    fprintf('\n')
    fprintf('Parameter schlecht gewählt !\n')
    fprintf('Rakete kommt nicht auf einen FPA von 0.\n')
    fprintf('Kreisbahn kann nur durch Kurskorrektur und Abbremsung erreicht werden.\n')
    fprintf('\n')
    PS = 3;
else
    fprintf('\n')
    fprintf('Parameter Ok - Kreisbahn kann erreicht werden.\n')
    fprintf('\n')
    PS = 5;
end

% Greenwich Nullmeridian beim Start  für Markierung der Abbildung
latG = deg2rad(linspace(-90,90,181)); 
lonG = 0;
xG0 = RE*cos(latG)*cos(lonG)/1000;
yG0 = RE*cos(latG)*sin(lonG)/1000;
zG0 = RE*sin(latG)/1000;

% Äquator
lonA = deg2rad(linspace(-180,180,361)); 
xA = RE*cos(lonA)/1000;
yA = RE*sin(lonA)/1000;
zA = 0*sin(lonA)/1000;


if PS ==5
    % Auf Kreisbahn einschwenken (auf ausreichende Geschwindigkeit kommen)
    AB4       = [vx4(end);vy4(end);vz4(end);x4(end);y4(end);z4(end);m4(end)];
    tspan4    = linspace(t4(end),t4(end)+tBA(3)*fracBurn,501);
    opts      = odeset('AbsTol',1.e-09,'RelTol',1.e-08,'Events', ...
        @(t5,Y)vkreis(t5,Y,P1));
    [t5,Y,te,ye,ie] = ode45(@(t5,Y) DGL_Thrust3Dtangential(t5,Y,P1), tspan4, ...
        AB4, opts);
    vx5(:)    = Y(:,1);
    vy5(:)    = Y(:,2);
    vz5(:)    = Y(:,3);
    x5(:)     = Y(:,4);
    y5(:)     = Y(:,5);
    z5(:)     = Y(:,6);
    m5(:)     = Y(:,7);
    h5(:)     = vecnorm([x5;y5;z5])-RE;
    
    % Freie Bewegung
    tUmlauf   = 2*pi/sqrt(muE/(RE+h5(end))^3);
    AB5       = [vx5(end);vy5(end);vz5(end);x5(end);y5(end);z5(end);m5(end)];
    tspan5    = linspace(t5(end),tUmlauf,501);
    opts      = odeset('AbsTol',1.e-09,'RelTol',1.e-08);
    [t6,Y]    = ode45(@(t6,Y) DGL_FreeTraj3D(t6,Y,P1), tspan5, AB5, opts);
    vx6(:)    = Y(:,1);
    vy6(:)    = Y(:,2);
    vz6(:)    = Y(:,3);
    x6(:)     = Y(:,4);
    y6(:)     = Y(:,5);
    z6(:)     = Y(:,6);
    m6(:)     = Y(:,7);
    h6(:)     = vecnorm([x6;y6;z6])-RE;

    % Greenwich Nullmeridian bei t=tUmlauf für Markierung der Abbildung
    lonG = OmegaE*tUmlauf;
    xGE = RE*cos(latG)*cos(lonG)/1000;
    yGE = RE*cos(latG)*sin(lonG)/1000;
    zGE = zG0;
    
    % Startpunkt bei t=tUmlauf für Markierung der Abbildung
    xSE = RE*cos(lat)*cos(lon+OmegaE*tUmlauf);
    ySE = RE*cos(lat)*sin(lon+OmegaE*tUmlauf);
    zSE = RE*sin(lat);

end


%%--------------------------------------------------------------------------
% Graphische Ausgabe
%--------------------------------------------------------------------------

figure('name','3D-Trajectory-Plot Ariane-1')
[XE,YE,ZE]=sphere(30);
surf(XE*RE/1000,YE*RE/1000,ZE*RE/1000,'linestyle',':')
alpha 0.75
LW =5;
hold on
plot3(xG0,yG0,zG0,'linewidth',2,'color',Colors(1,:),'linestyle',':');
plot3(xA,yA,zA,'linewidth',2,'color',Colors(1,:),'linestyle',':');
plot3(x1/1000,y1/1000,z1/1000,'color',Colors(2,:),'LineWidth',LW);
plot3(x2/1000,y2/1000,z2/1000,'color',Colors(4,:),'LineWidth',LW);
plot3(x3/1000,y3/1000,z3/1000,'color',Colors(3,:),'LineWidth',LW);
plot3(x4/1000,y4/1000,z4/1000,'color',Colors(6,:),'LineWidth',LW);
if PS == 5
    plot3(xGE,yGE,zGE,'linewidth',2,'color',Colors(1,:),'linestyle',':');
    plot3(x5/1000,y5/1000,z5/1000,'color',Colors(7,:),'LineWidth',LW);
    plot3(x6/1000,y6/1000,z6/1000,'color',Colors(9,:),'LineWidth',LW);
end
% Abschussort bei Zeitpunkt t=0
plot3(x0/1000,y0/1000,z0/1000,   '+y','markersize' ,6,'linewidth',2,...
      'Color',Colors(10,:));
%Greenwhich zum Zeitpunkt Start
plot3(RE*cosd(50)/1000,0,RE*sind(50)/1000,'+y','markersize',6,'linewidth',2) 
if PS == 5 
    % Abschussort bei Zeitpunkt t=tUmlauf
    plot3(xSE/1000,ySE/1000,zSE/1000,'+y','markersize' ,6,'linewidth',2,...
          'Color',Colors(10,:));
    %Greenwhich zum Zeitpunkt  t=tUmlauf
    lonG = OmegaE*tUmlauf;
    xGEP = RE*cosd(50)*cos(lonG)/1000;
    yGEP = RE*cosd(50)*sin(lonG)/1000;
    zGEP = RE*sind(50)/1000;
    plot3(xGEP,yGEP,zGEP,'+y','markersize' ,6,'linewidth',2) 
end

axis equal 
grid on
str= "3D-Trajectory-Plot Ariane-1";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
xlabel('x in km','FontSize',14); 
zlabel('z in km','FontSize',14); 
ylabel('y in km','FontSize',14)
set(gca,'FontSize',16);


%---------------------------------------------------------------------
% Übersichts-Plots Startphase
taxend = sum(tBA)+60;  %Zeitdauer für den Plot
taxend = 900;

figure('name','Übersicht-Plots ohne freie Trajektorie')
subplot(2,2,1)
hold on
p(1)=plot(t1,h1/1000,'color',Colors(2,:),'Linewidth',2);
p(2)=plot(t2,h2/1000,'color',Colors(4,:),'Linewidth',2);
p(3)=plot(t3,h3/1000,'color',Colors(3,:),'Linewidth',2);
p(4)=plot(t4,h4/1000,'color',Colors(6,:),'Linewidth',2);
if PS== 5 
    p(5)=plot(t5,h5/1000,'color',Colors(7,:),'Linewidth',2);
    p(6)=plot(t6,h6/1000,'color',Colors(9,:),'Linewidth',2);
end
lgdstr(1,:)='Stufe 1 bis Pitch    ';
lgdstr(2,:)='Stufe 1 bis Burn-Out ';
lgdstr(3,:)='Stufe 2 bis Burn-Out ';
lgdstr(4,:)='Stufe 3 bis FPA=0    ';
lgdstr(5,:)='Stufe 3 bis Kreisbahn';
lgdstr(6,:)='Freie Trajektorie    ';
str= "h(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:PS+1),lgdstr(1:PS+1,:),'location','southeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, taxend, 0 h4(end)*1.2/1000]);
legend box off;
xlabel('t in s','FontSize',14); 
ylabel('h in km','FontSize',14)
grid on
set(gca,'FontSize',16);

subplot(2,2,2)
v1 = vecnorm([vx1;vy1;vz1]);
v2 = vecnorm([vx2;vy2;vz2]);
v3 = vecnorm([vx3;vy3;vz3]);
v4 = vecnorm([vx4;vy4;vz4]);
hold on
p(1)=plot(t1,v1/1000,'color',Colors(2,:),'Linewidth',2);
p(2)=plot(t2,v2/1000,'color',Colors(4,:),'Linewidth',2);
p(3)=plot(t3,v3/1000,'color',Colors(3,:),'Linewidth',2);
p(4)=plot(t4,v4/1000,'color',Colors(6,:),'Linewidth',2);
if PS== 5 
    v5 = vecnorm([vx5;vy5;vz5]);
    v6 = vecnorm([vx6;vy6;vz6]);
    p(5)=plot(t5,v5/1000,'color',Colors(7,:),'Linewidth',2);
    p(6)=plot(t6,v6/1000,'color',Colors(9,:),'Linewidth',2);
end
str= "v(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:PS+1),lgdstr(1:PS+1,:),'location','southeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, taxend, 0 , max(8,1.2*v4(end)/1000)]);
xlabel('t in s','FontSize',14); 
ylabel('v in km/s','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);

% Projektion der Geschwindigkeit (inklusive mitgenommene Erdrotation)
% auf die Tangentialebene unter der Rakete
% Der so berechnete Winkel Gamma ist somit der FPA im Inertialsystem unter 
% Berücksichtigung der Anfangsgeschwindigkeit durch die Erdrotation

vloc   = sqrt(vx1.^2 + vy1.^2 + vz1.^2);
Rloc   = sqrt(x1.^2 + y1.^2 + z1.^2);
gamma1 = rad2deg(pi/2-acos((x1.*vx1+y1.*vy1+z1.*vz1)./vloc./Rloc));
%
vloc   = sqrt(vx2.^2 + vy2.^2 + vz2.^2);
Rloc   = sqrt(x2.^2 + y2.^2 + z2.^2);
gamma2 = rad2deg(pi/2-acos((x2.*vx2+y2.*vy2+z2.*vz2)./vloc./Rloc));
%
vloc   = sqrt(vx3.^2 + vy3.^2 + vz3.^2);
Rloc   = sqrt(x3.^2 + y3.^2 + z3.^2);
gamma3 = rad2deg(pi/2-acos((x3.*vx3+y3.*vy3+z3.*vz3)./vloc./Rloc));
%
vloc   = sqrt(vx4.^2 + vy4.^2 + vz4.^2);
Rloc   = sqrt(x4.^2 + y4.^2 + z4.^2);
gamma4 = rad2deg(pi/2-acos((x4.*vx4+y4.*vy4+z4.*vz4)./vloc./Rloc));
%
if PS== 5 
    vloc   = sqrt(vx5.^2 + vy5.^2 + vz5.^2);
    Rloc   = sqrt(x5.^2 + y5.^2 + z5.^2);
    gamma5 = rad2deg(pi/2-acos((x5.*vx5+y5.*vy5+z5.*vz5)./vloc./Rloc));
    vloc   = sqrt(vx6.^2 + vy6.^2 + vz6.^2);
    Rloc   = sqrt(x6.^2 + y6.^2 + z6.^2);
    gamma6 = rad2deg(pi/2-acos((x6.*vx6+y6.*vy6+z6.*vz6)./vloc./Rloc));
end

subplot(2,2,3)
hold on
p(1)=plot(t1,gamma1,'color',Colors(2,:),'Linewidth',2);
p(2)=plot(t2,gamma2,'color',Colors(4,:),'Linewidth',2);
p(3)=plot(t3,gamma3,'color',Colors(3,:),'Linewidth',2);
p(4)=plot(t4,gamma4,'color',Colors(6,:),'Linewidth',2);
if PS==5 
   p(5)=plot(t5,gamma5,'color',Colors(7,:),'Linewidth',2);
   p(6)=plot(t6,gamma6,'color',Colors(9,:),'Linewidth',2);
end
hp2=legend(p(1:PS+1),lgdstr(1:PS+1,:),'location','northeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
str= "\gamma(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
axis([0, taxend, -10 40]);
xlabel('t in s','FontSize',14); 
ylabel('\gamma','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);

subplot(2,2,4)
semilogy(t1,m1,'color',Colors(2,:),'Linewidth',2);
hold on
p(1)=semilogy(t1,m1,'color',Colors(2,:),'Linewidth',2);
p(2)=semilogy(t2,m2,'color',Colors(4,:),'Linewidth',2);
p(3)=semilogy(t3,m3,'color',Colors(3,:),'Linewidth',2);
p(4)=semilogy(t4,m4,'color',Colors(6,:),'Linewidth',2);
if PS==5 
    p(5)=semilogy(t5,m5,'color',Colors(7,:),'Linewidth',2);
    p(6)=semilogy(t6,m6,'color',Colors(9,:),'Linewidth',2);
end
tall = linspace(0,taxend);
semilogy(tall,mLA*tall./tall,'color',Colors(8,:),'Linewidth',2);
htext=text(10,1.25*mLA,'Nutzlast','color',Colors(8,:));
set(htext,'FontSize',14,'FontWeight','normal'); 
str= "m(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:PS+1),lgdstr(1:PS+1,:),'location','northeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, taxend, 1e3 2*m0]);
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
St     = P1.Stage;
R      = norm([Y(4) Y(5) Y(6)]);
rho    = P1.rho0*exp(-(R-P1.RE)/P1.H0);
% Geschwindigkeit relativ zur Atmosphäre für Reibung
vE_rot = P1.OmegaE*sqrt(Y(4)^2+Y(5)^2);
vRv    = [Y(1)+vE_rot*sin(atan2(Y(5),Y(4))) ...
         Y(2)-vE_rot*cos(atan2(Y(5),Y(4))) Y(3)];
vR     = norm(vRv);
F_D    = 0.5*rho*P1.cwD*vR^2;
C1     = P1.FSA(St)/Y(7);
C2 = -F_D/vR/Y(7);
C3     = -P1.muE/R^3;
dY     = [C1*P1.S(1)+C2*vRv(1)+C3*Y(4);        
          C1*P1.S(2)+C2*vRv(2)+C3*Y(5);  
          C1*P1.S(3)+C2*vRv(3)+C3*Y(6);  
          Y(1);
          Y(2);
          Y(3);
          -P1.dmA(St)];
end

function dY = DGL_Thrust3D(t, Y, P1)
% Diese Variante der Differentialgleichung geht davon aus, dass
% die Rakete sich entlang ihrer Geschwindigkeit ausrichtet, wobei
% die Geschwindigkeit, die sie aufgrund der Erdrotation mitnimmt
% nicht berücksichtigt wird.
% Y(1):vx, Y(2):vy, Y(3):vz, Y(4):x Y(5):y Y(6):z Y(7):m
St     = P1.Stage;
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
C1     = P1.FSA(St)/vRS/Y(7);
C2 = -F_D/vR/Y(7);
C3     = -P1.muE/R^3;
dY     = [C1*vRSv(1)+C2*vRv(1)+C3*Y(4);        
          C1*vRSv(2)+C2*vRv(2)+C3*Y(5);  
          C1*vRSv(3)+C2*vRv(3)+C3*Y(6);  
          Y(1);
          Y(2);
          Y(3);
          -P1.dmA(St)];
end

function dY = DGL_Thrust3Dtangential(t, Y, P1)
% Diese Variante der Differentialgleichung sorgt für eine Beschleunigung
% tangential zur Erdoberfläche, Projektion der Geschwindigkeit auf
% die Tangentialebene zum Erdboden.
% Y(1):vx, Y(2):vy, Y(3):vz, Y(4):x Y(5):y Y(6):z Y(7):m
St     = P1.Stage;
R      = norm([Y(4) Y(5) Y(6)]);
rho    = P1.rho0*exp(-(R-P1.RE)/P1.H0);
% Geschwindigkeit relativ zur Atmosphäre für Reibung
vE_rot = P1.OmegaE*sqrt(Y(4)^2+Y(5)^2);
vRv    = [Y(1)+vE_rot*sin(atan2(Y(5),Y(4))) ...
    Y(2)-vE_rot*cos(atan2(Y(5),Y(4))) Y(3)];
vR     = norm(vRv);
% Geschwindigkeit tangential zur Erdoberfläche
vSenk  = (Y(1)*Y(4)+Y(2)*Y(5)+Y(3)*Y(6))/R;
vTanv  = [Y(1)-vSenk*Y(4) Y(2)-vSenk*Y(5) Y(3)-vSenk*Y(6)]/R;
vTan   = norm(vTanv);
F_D    = 0.5*rho*P1.cwD*vR^2;
C1     = P1.FSA(St)/vTan/Y(7);
C2 = -F_D/vR/Y(7);
C3     = -P1.muE/R^3;
dY     = [C1*vTanv(1)+C2*vRv(1)+C3*Y(4);        
          C1*vTanv(2)+C2*vRv(2)+C3*Y(5);  
          C1*vTanv(3)+C2*vRv(3)+C3*Y(6);  
          Y(1);
          Y(2);
          Y(3);
          -P1.dmA(St)];
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

% Halte an, wenn eine Geschwindigkeit erreicht wird, die die gewünschte
% Kreisbahn ergibt.
function [value,isterminal,direction] = vkreis(t, Y, P1)
        vloc2 = Y(1)^2 + Y(2)^2 + Y(3)^2;
        Rloc  = sqrt(Y(4)^2 + Y(5)^2 + Y(6)^2);
        value = vloc2-P1.muE/Rloc;
                            % identifiziere passende Geschwindkgkeit
        isterminal = 1;     % Integration anhalten
        direction = 0;      % Geschwindigkeit = 0 wird von unten erreicht
end
