% -------------------------------------------------------------------------
% IonenantriebMond.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet numerisch (ODE45) die Bahnen/Flugzeiten
% für einen Ionenantrieb aus einem LEO zum Mond 
% im Vergleich zum Hohmannn-Transfer 
%
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


% Parameter
GME = 398600.4415;   % G*ME in km^3/s^2
MM  = 7.348e22;      % Mondmasse in kg
ME  = 5.972e24;      % Erdmasse in kg
GMM = GME*MM/ME;     % G*MM in km^3/s^2
RE  = 6378;          % Erdradius in km
rEM = 384000;        % Abstand Erde Mond
rEM0= rEM-58000;     % Lagrange Punkt zw. Erde und Mond
rM  = 6000;          % Radius Mondbahn
r2  = rEM+rM;        % Mondbahn + Höhe Mondumlaufbahn
r20 = r2;

% Hilfsgrößen
phi = linspace(0,2*pi,361);
xK    = cos(phi);
yK    = sin(phi);

% Raumsonde Parameter
h1  = 1000;          % Höhe LEO Bahn
r1  = h1+RE;         % geozentrischer Radius Ausgangsbahn
r10 = r1;            % geozentrischer Radius Ausgangsbahn
v10 = sqrt(GME/r10); % Anfangsgeschwindigkeit                 
m0 = 10000;          % Anfangsmasse in kg
FS = 0.05;           % Anfangsschub Ionentriebwerk in N

% Parameter für ODE45
P1.m0   = m0;
P1.GME  = GME;
P1.GMM  = GME*MM/ME;
P1.FS   = FS;
P1.mend = 2*m0/3;
P1.rEM  = rEM;

%% Berechnung des Geschwindigkeitszuwaches und der Bahn nach Hohmann

aH0 = 0.5*(r10+r20);           % Transferellipse Halbachse
vP = sqrt(GME*r2/aH0/r1);      % v bei Perigäum
vA = sqrt(GME*r1/aH0/r20);     % v bei Apogäum
Delv  = vP-v10;                % Geschwindigkeitsbudget in km/s
tH0= pi*sqrt(aH0^3/GME)/3600;  % Flugzeit für Hohmann Ellipse in Stunden
eccH0 = sqrt(1-r1*r2/aH0^2);   % Exzentrizität

% Betriebszeit des Ionentriebwerks auf Basis Hohmann-Transferzeit
tBurn = 3.75*tH0*3600;
P1.stopBurn  = tBurn;

%% Bahnform mit Zeitberechnung über Lösung Keplergleichung
% u = linspace(0,180,3601);
% % cosd, sind = Cosinus, Sinus in Grad (degree)
% cosu=cosd(u);
% sinu=sind(u);
% gamma = aH0*(1-eccH0.*eccH0);
% % Bahnkoordinaten 
% xH = gamma*(cosu./(1+eccH0.*cosu));
% yH = gamma*(sinu./(1+eccH0.*cosu));  
% rH = sqrt(xH.^2+yH.^2);

% %% Numerische Berechnung Hohmann-Ellipse über Zeit
% 
% %Parameter für Routine Bahnberechnung
% tH = linspace(0,tH0*3600,1000);
% BaPaK.eP  = eccH0;  %Exzentrizität
% BaPaK.T0  = 0;      %Perihelzeit
% BaPaK.aP  = aH0;    %große Halbachse
% BaPaK.qp  = aH0*(1-eccH0);  %Perihelabstand
% BData=BahnBerechung(GME, tH, BaPaK);  %siehe Kapitel Himmelsmechanik
% xH1 = BData.xyz(1,:);
% yH1 = BData.xyz(2,:);
% rH1 = sqrt(xH1.^2+yH1.^2);

%% Numerische Berechnung Ionenantrieb Abstand Erde 

tend = 5*tH0*3600;
t  = linspace(0,tend,10000);
x0  = RE+h1;
y0  = 0; 
vx0 = 0; 
vy0 = v10;
AB = [x0;vx0;y0;vy0];          %AB für DGL Ionentriebwerk
% MATLAB's Runge-Kutta ode45 Routine Ionentriebwerk
opts = odeset('AbsTol',1.e-9,'RelTol',1.e-8,'events',@MyEvent1);
[tS,YS, TES, YE, IE]=ode45(@DGL_IonTriebwerk, t, AB,opts,P1);

xS = YS(:,1);
yS = YS(:,3);
rS = sqrt(xS.^2+yS.^2);

%Berechnung Ionentriebwerksschub und Masse Raumsonde als Funktion der Zeit
m     = max(P1.mend,P1.m0 - (P1.m0  - P1.mend).*tS/tBurn);
FSist = P1.FS*(1-tS.^2/tBurn.^2);
FSist = max(0,FSist);  
ratioF = FSist./(GME./rS.^2);  %Verhältnis zur Gravitationskraft

%% Numerische Berechnung Hohmann-Transfer als Funktion der Zeit

vy0 = sqrt(GME/x0)+ Delv;      %Impulsschub mit Delv @ t=0          
AB = [x0;vx0;y0;vy0];          %AB für DGL Hohmann
t  = linspace(0,1.2*tH0*3600,10000);

% MATLAB's Runge-Kutta ode45 Routine Hohmann-like-Ellipse
opts = odeset('AbsTol',1.e-9,'RelTol',1.e-8,'events',@MyEvent1);
[tH,YH, TEH, YE, IE]=ode45(@DGL_Hohmann_like, t, AB,opts,P1);

xH = YH(:,1);
yH = YH(:,3);
rH = sqrt(xH.^2+yH.^2);

%%  Graphische Darstellung

figure('Name','Abstand von Erde und Schub über Zeit für Ionenantrieb')
yyaxis left
h(1)=plot(tS/3600, rS/1000,'color',Colors(3,:),'Linewidth',2,...
         'LineStyle',Style(1));
hold on
h(2)=plot(tH/3600, rH/1000,'color',Colors(2,:),'Linewidth',2,...
         'LineStyle',Style(1));
h(3)=line([0 TES/3600], [rEM/1000 rEM/1000],...
      'color',Colors(10,:),'Linewidth',2,'LineStyle',Style(1));
h(4)=line([0 TES/3600], [rEM0/1000 rEM0/1000],...
      'color',Colors(10,:),'Linewidth',2,'LineStyle',Style(3));
h(5)=line([TES TES]/3600, [0 rEM/1000],...
     'color',Colors(3,:),'Linewidth',1,'LineStyle',Style(3));
h(6)=line([TEH TEH]/3600, [0 rEM/1000],...
     'color',Colors(2,:),'Linewidth',1,'LineStyle',Style(3));
grid on
ylabel('r in 1000 km')
xlim([0 TES*1.1]/3600); ylim([RE/1000 1.1*r20/1000]);
yyaxis right
h(7)= plot(tS/3600, FSist,'color',Colors(5,:),...
     'Linewidth',1,'LineStyle',Style(1));
hold on
xlabel('t in h')
ylabel('F_S in N ')
ylim([0 max(FSist)*1.1]);
legend(h,' Ionenantrieb',' Hohmann-Bahn',' Mondbahn',...
         ' Lagrange-Pkt. L1',' t_{ Ionenantrieb}',...
         ' t_{ Hohmann}',' Schub','location','eastoutside');   
legend box off
ttl=title('Ionenantrieb vs Hohmann-Transfer r(t)');
set(gca,'FontSize',14);
set(ttl,'FontSize',14,'FontWeight','normal');

%%

figure('Name','Zeit über Abstand Erde und Gravitation/Schub ')
yyaxis left
h(1)=plot(rS/1000,tS/3600,'color',Colors(3,:),'Linewidth',2,...
         'LineStyle',Style(1));
hold on
h(2)=plot(rH/1000, tH/3600, 'color',Colors(2,:),'Linewidth',2,...
         'LineStyle',Style(1));
h(3)=line([rEM/1000 rEM/1000], [0 TES/3600],...
       'color',Colors(10,:),'Linewidth',1,'LineStyle',Style(1));
h(4)=line([rEM0/1000 rEM0/1000],[0 TES/3600], ...
       'color',Colors(10,:),'Linewidth',1,'LineStyle',Style(3));
h(5)=line([0 rEM/1000],[TES TES]/3600, ...
     'color',Colors(3,:),'Linewidth',1,'LineStyle',Style(3));
h(6)=line([0 rEM/1000],[TEH TEH]/3600,...
     'color',Colors(2,:),'Linewidth',1,'LineStyle',Style(3));
grid on
xlabel('r in 1000 km')
ylabel('t in h')
ylim([0 TES*1.1]/3600); xlim([RE/1000 1.1*r20/1000]);
yyaxis right
h(7)= plot(rS/1000, FSist,'color',Colors(5,:),...
     'Linewidth',1,'LineStyle',Style(1));
hold on
xlabel('r in 1000 km')
ylabel('F_S in N ')
ylim([0 max(FSist)*1.1]);
legend(h,' Ionenantrieb',' Hohmann-Bahn',' Mondbahn',...
         ' Lagrange-Pkt. L1',' t_{ Ionenantrieb}',...
         ' t_{ Hohmann}',' Schub','location','eastoutside');   
legend box off
ttl=title('Ionenantrieb vs Hohmann-Transfer t(r)');
set(gca,'FontSize',14);
set(ttl,'FontSize',14,'FontWeight','normal');


%% Trajektorien-Graphik

figure('Name','Trajektorien')
plot(xS/1000,yS/1000,'color',Colors(3,:),...
     'Linewidth',2,'LineStyle',Style(3));
hold on
plot(xH/1000,yH/1000,'color',Colors(2,:),...
     'Linewidth',2,'LineStyle',Style(2));
plot(r10*xK/1000,r10*yK/1000,'color',Colors(3,:),...
     'Linewidth',1,'LineStyle',Style(1));
plot((-rEM+rM*xK)/1000,rM*yK/1000,'color',Colors(10,:),...
     'Linewidth',1,'LineStyle',Style(1));
plot(xK*rEM/1000, yK*rEM/1000,'color',Colors(10,:),...
     'Linewidth',1,'LineStyle',Style(1));
plot(xK*rEM0/1000, yK*rEM0/1000,'color',Colors(10,:),...
     'Linewidth',1,'LineStyle',Style(2));
grid on
axis equal
xlim([-rEM rEM]*1.1e-3); ylim([-rEM rEM]*1.1e-3);
legend(' Ionenantrieb',' Hohmann-Bahn',' LEO', ' Mondbahn', ' L1 ');   
legend box off
ttl=title('Trajektorien');
set(ttl,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14);



%% ------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen
 
% DGL
function dY = DGL_IonTriebwerk(t, Y, P1)
% Y(1):x, Y(2):vx, Y(3):y, Y(4):vy
rEM = P1.rEM;
tBurn = P1.stopBurn;
m     = max(P1.mend,P1.m0 - (P1.m0  - P1.mend).*t/tBurn);
FSist = P1.FS*(1-t.^2/tBurn.^2);
as = max(0,FSist)./m;  
r   = sqrt(Y(1).^2+Y(3).^2);
rM  = sqrt((Y(1)-rEM).^2+Y(3).^2);
phi = atan2(Y(3),Y(1));
dY = [Y(2);
      -P1.GME*Y(1)./r^3  - P1.GMM*(Y(1)-rEM)./rM^3 - as*sin(phi); 
      Y(4);
      -P1.GME*Y(3)./r^3 - P1.GMM*Y(3)./rM^3 + as*cos(phi)];  
end 

% DGL
function dY = DGL_Hohmann_like(t, Y, P1)
% Y(1):x, Y(2):vx, Y(3):y, Y(4):vy
rEM = P1.rEM;
r   = sqrt(Y(1).^2+Y(3).^2);
rM  = sqrt((Y(1)-rEM).^2+Y(3).^2);
dY = [Y(2);
      -P1.GME*Y(1)./r^3  - P1.GMM*(Y(1)-rEM)./rM^3; 
      Y(4);
      -P1.GME*Y(3)./r^3 - P1.GMM*Y(3)./rM^3];  
end 

function BData=BahnBerechung(GM, t, BaPaK)
    %   qPn      Periheldistanz 
    %   aPn      grioße Halbachse 
    %   ePn      ExzentrizitÃ¤t
    ePn=BaPaK.eP;
    aPn=BaPaK.aP;
    qPn=BaPaK.qp;
    tau=sqrt(GM)*t;

    % Bahnkoordinaten aus Stumpff-Funktionen Näherung fÃ¼r Parabelbahn
    fac =0.5*ePn;
    E2=100;
    E20=0;
    iter=0;
    kf=sqrt(GM/(qPn*(1+ePn)));
    while abs(E2-E20)>eps
      iter=iter+1;
      E20 = E2;
      A = 1.5*sqrt(fac/(qPn*qPn*qPn)).*tau;  
      B = (sqrt(A.*A+1.0)+A).^(1/3);
      u  = B - 1.0./B;  
      u2 = u.*u;  
      E2 = u2*(1.0-ePn)./fac;
      [c1,c2,c3]=Stumpff(E2); 
      fac = 3.0*ePn.*c3;
    end
    R  = qPn*(1 + u2.*c2*ePn./fac );
    x=qPn*(1.0-u2.*c2./fac);
    y=qPn*sqrt((1.0+ePn)./fac).*u.*c1;
    z=0*x;
    vx=-kf*y./R;
    vy=kf*(x./R+ePn); 
    BData.xyz = [x;y;z]; 
    BData.v   = [vx; vy; 0.0*vx];
    BData.Time = t;
end

function [value,isterminal,direction] = MyEvent1(t,YS,P1)
    %Ereignisfunktion bis Eintreten in Mondumlaufbahn
    value      = P1.rEM - sqrt(YS(1).^2+YS(3).^2); % detect distance
    isterminal = 1;        % stop the integration
    direction  = 0;        % negative direction
end

%Ende Funktionen
% -------------------------------------------------------------------------

