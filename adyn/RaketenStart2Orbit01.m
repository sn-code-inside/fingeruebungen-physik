% -------------------------------------------------------------------------
% RakentenStart2Orbit01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet Start einer dreistufigen Rakete von der Erdoberfläche
% unter Berücksichtigung von Reibung und Nickwinkel-Korrekturen.
% Berechnung im mitbewegten KOS.
% 
% Beispiele nach Messerschmidt "Raumfahrtsystem, Springer Verlag 2017
% S. 141
% Alle Rechnungen in kg m s
% Der Bogen s(t) wird durch x(t) angenähert.
%--------------------------------------------------------------------------


%%
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;

%%
% Initialisierung
% alle Daten in kg, m und s 
G   = 6.67430e-11;                 % G in m^3 /kg /s2
ME  = 5.98e24;                     % Masse der Erde in kg
RE  = 6.378e6;                     % Erdradius in m
gE  = G*ME/RE^2;                   % gE in m/s^2
muE = G*ME;

% Parameter der Raketen

% Daten Rakete A 

mSA     = [87.0 10.0 3.0]*1e3;
mTA     = [80.0	8.0 1.320]*1e3;
m0      = sum(mSA);
m0A     = [m0	m0-mSA(1)	m0-mSA(1)-mSA(2)];
tBA     = [160	80	66];
ISpezA  = [300 300 300];
dmA     = mTA./tBA;
FSA     = ISpezA.*dmA*gE;
mLA     = mSA(3)- mTA(3);
sigmaL  = mLA/m0;


% Parameter der DGL
rho0       = 1.752;               % Dichte in kg/m^3
H0         = 6.7e3;               % H0 in m Barometr. Höhenformel
cwD        = 1.5;                 % Widerstandsbeiwert x Fläche in m^2 
hT         = 120e3;               % Übergangshöhe in m
hK         = 200e3;               % Höhe finale Kreisbahn in m
tCoast     = 25;                   % Dauer Coasting Phase in s

P1.gE       = gE;                 
P1.muE      = muE;
P1.H0       = H0;                    
P1.rho0     = rho0;
P1.RE       = RE;
P1.cwD      = cwD;
P1.dmA      = dmA;
P1.tBA      = tBA;
P1.FSA      = FSA;

% Anfangswerte 
v0     = 0.001;
gamma0 = pi/2;
h0     = 0;
x0     = 0;
m0     = m0A(1);

% Allgemeiner Option-Set für MATLABs Runge-Kutta ode45 Routine 
opts   = odeset('AbsTol',1.e-09,'RelTol',1.e-08);

% Berechnung Anfangsphase bis Nickkorrektur in Höhe von 1km 
optsNick  = odeset('AbsTol',1.e-9,'RelTol',1.e-8,'Events',@myeventNick);
AB0       = [v0;gamma0;h0;x0;m0];   % AB für DGL
tNick     = 20;                     % Zeit bis späteste Nickkorrektur in s
P1.Stage  = 1;
tspan1    = linspace(0,tNick,101);
[t1,Y,te,ye,ie] = ode45(@(t1,Y) DGL_Thrust(t1,Y,P1), tspan1, AB0,optsNick);
v1(:)     = Y(:,1);
gamma1(:) = Y(:,2);
h1(:)     = Y(:,3);
x1(:)     = Y(:,4);
m1(:)     = Y(:,5);

% Berechnung erste Stufe nach Nickkorrektur bis Burn-Out Stufe 1
Del_gamma_v= deg2rad([0.5,1,5.5,6.5]); % Auswahl möglicher Nickkorrekturen
Del_gamma  = Del_gamma_v(4);
AB1       = [v1(end);gamma1(end)-Del_gamma;h1(end);x1(end);m1(end)];    
P1.Stage  = 1;
if te < tNick 
    tNick = te;
end
tspan1    = linspace(tNick,tBA(1),501);
[t2,Y]    = ode45(@DGL_Thrust,tspan1, AB1,opts,P1);
v2(:)     = Y(:,1);
gamma2(:) = Y(:,2);
h2(:)     = Y(:,3);
x2(:)     = Y(:,4);
m2(:)     = Y(:,5);

% Berechnung zweite Stufe bis Burn-Out Stufe 2 und evetl caosting phase 
P1.Stage  = 2;
AB2       = [v2(end);gamma2(end);h2(end);x2(end);m0A(2)];   
tspan2    = linspace(t2(end),t2(end)+tBA(2),501);
[t3a,Y]    = ode45(@DGL_Thrust,tspan2, AB2,opts,P1);
v3a(:)     = Y(:,1);
gamma3a(:) = Y(:,2);
h3a(:)     = Y(:,3);
x3a(:)     = Y(:,4);
m3a(:)     = Y(:,5);
if tCoast >0
%Coasting
    AB3       = [v3a(end);gamma3a(end);h3a(end);x3a(end);m0A(3)];   
    tspan3    = linspace(t3a(end),t3a(end)+tCoast,501);
    [tc,Y]    = ode45(@DGL_FreeTraj, tspan3, AB3, opts, P1);
    t3        = [t3a; tc];
    v3(:)     = [v3a transpose(Y(:,1))];
    gamma3(:) = [gamma3a transpose(Y(:,2))];
    h3(:)     = [h3a Y(:,3).'];
    x3(:)     = [x3a Y(:,4).'];
    m3(:)     = [m3a Y(:,5).'];
else
    t3        = t3a;
    v3(:)     = v3a;
    gamma3(:) = gamma3a;
    h3(:)     = h3a;
    x3(:)     = x3a;
    m3(:)     = m3a;
end 

% Berechnung dritte Stufe bis Burn-Out Stufe 3
AB3       = [v3(end);gamma3(end);h3(end);x3(end);m0A(3)];    
% opts3     = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'Events',@myeventGamma);
P1.Stage  = 3;
fracBurn  = 1.0;    %Anteil des Treibstoffs von Stufe 
tspan3    = linspace(t3(end),t3(end)+tBA(3)*fracBurn,501);
% [t3,Y,te,ye,ie] = ode45(@(t3,Y) DGL_Launch(t3,Y,P1), tspan3, AB3, opts3);
[t4,Y]    = ode45(@DGL_Thrust, tspan3, AB3,opts,P1);
v4(:)     = Y(:,1);
gamma4(:) = Y(:,2);
h4(:)     = Y(:,3);
x4(:)     = Y(:,4);
m4(:)     = Y(:,5);

%Phase nach Burn-Out, free trajectory
AB4       = [v4(end);gamma4(end);h4(end);x4(end);m4(end)];   
opts4     = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'Events',@myeventxend);
tspan4    = linspace(t4(end),t4(end)+300,501);
[t5,Y,te,ye,ie] = ode45(@(t5,Y) DGL_FreeTraj(t5,Y,P1), tspan4, AB4, opts4);
v5(:)     = Y(:,1);
gamma5(:) = Y(:,2);
x5(:)     = Y(:,4);
h5(:)     = Y(:,3);
m5(:)     = Y(:,5);

%% Graphische Ausgabe
%--------------------------------------------------------------------------

figure('name','h-Plots')
hold on
p(1)=plot(x1/1000,h1/1000,'color',Colors(2,:),'Linewidth',2);
p(2)=plot(x2/1000,h2/1000,'color',Colors(8,:),'Linewidth',2);
p(3)=plot(x3/1000,h3/1000,'color',Colors(4,:),'Linewidth',2);
p(4)=plot(x4/1000,h4/1000,'color',Colors(3,:),'Linewidth',2);
p(5)=plot(x5/1000,h5/1000,'color',Colors(7,:),'Linewidth',2);
for k=1:5 
    lgdstr(k,:) = strcat("Stufe ",num2str(k-1));
end    
lgdstr(1,:)='Stufe 1 bis Pitch';
lgdstr(5,:)='Freie Trajektorie';
grid on
str= "h(x)-Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'location','southeast');  
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, inf, 0 inf]);
xlabel('x in km','FontSize',14); 
ylabel('h in km','FontSize',14)
legend box off;
set(gca,'FontSize',16);

figure('name','Übersicht-Plots')
subplot(2,2,1)
hold on
p(1)=plot(t1,h1/1000,'color',Colors(2,:),'Linewidth',2);
p(2)=plot(t2,h2/1000,'color',Colors(8,:),'Linewidth',2);
p(3)=plot(t3,h3/1000,'color',Colors(4,:),'Linewidth',2);
p(4)=plot(t4,h4/1000,'color',Colors(3,:),'Linewidth',2);
p(5)=plot(t5,h5/1000,'color',Colors(7,:),'Linewidth',2);
str= "h(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'location','southeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, t5(end), 0 inf]);
legend box off;
xlabel('t in s','FontSize',14); 
ylabel('h in km','FontSize',14)
grid on
set(gca,'FontSize',16);

subplot(2,2,2)
hold on
p(1)=plot(t1,v1/1000,'color',Colors(2,:),'Linewidth',2);
p(2)=plot(t2,v2/1000,'color',Colors(8,:),'Linewidth',2);
p(3)=plot(t3,v3/1000,'color',Colors(4,:),'Linewidth',2);
p(4)=plot(t4,v4/1000,'color',Colors(3,:),'Linewidth',2);
p(5)=plot(t5,v5/1000,'color',Colors(7,:),'Linewidth',2);
str= "v(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'location','southeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, t5(end), 0 inf]);
xlabel('t in s','FontSize',14); 
ylabel('v in km/s','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);

subplot(2,2,3)
hold on
p(1)=plot(t1,rad2deg(gamma1),'color',Colors(2,:),'Linewidth',2);
p(2)=plot(t2,rad2deg(gamma2),'color',Colors(8,:),'Linewidth',2);
p(3)=plot(t3,rad2deg(gamma3),'color',Colors(4,:),'Linewidth',2);
p(4)=plot(t4,rad2deg(gamma4),'color',Colors(3,:),'Linewidth',2);
p(5)=plot(t5,rad2deg(gamma5),'color',Colors(7,:),'Linewidth',2);
str= "\gamma(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:k),lgdstr(1:k,:),'location','northeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, t5(end), -30 90]);
xlabel('t in s','FontSize',14); 
ylabel('\gamma in Â°','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);


subplot(2,2,4)
p(1)=semilogy(t1,m1,'color',Colors(2,:),'Linewidth',2);
hold on
p(2)=semilogy(t2,m2,'color',Colors(8,:),'Linewidth',2);
p(3)=semilogy(t3,m3,'color',Colors(4,:),'Linewidth',2);
p(4)=semilogy(t4,m4,'color',Colors(3,:),'Linewidth',2);
p(5)=semilogy(t5,m5,'color',Colors(7,:),'Linewidth',2);
str= "m(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:k),lgdstr(1:k,:),'location','southwest'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, t5(end), 1e2 m0]);
xlabel('t in s','FontSize',14); 
ylabel('m in kg','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);


%% Funktionen
%--------------------------------------------------------------------------
function dY = DGL_Thrust(t, Y, P1)
% Y(1):v, Y(2):gamma, Y(3):h, Y(4):x Y(5): m
St    = P1.Stage;
rho   = P1.rho0*exp(-Y(3)/P1.H0); 
F_S   = P1.FSA(St);  
F_D   = 0.5*rho*P1.cwD*Y(1)^2;
g     = P1.muE/(P1.RE+Y(3))^2;
dY    = [(F_S-F_D)/Y(5) - g*sin(Y(2));
        -(g/Y(1)-Y(1)/(P1.RE+Y(3)))*cos(Y(2));
         Y(1)*sin(Y(2));
         Y(1)*cos(Y(2));
         -P1.dmA(St)];
end

function dY = DGL_FreeTraj(t, Y, P1)
% Y(1):v, Y(2):gamma, Y(3):h, Y(4):x Y(5): m
rho   = P1.rho0*exp(-Y(3)/P1.H0); 
F_S   = 0;  
F_D   = 0.5*rho*P1.cwD*Y(1)^2;
g     = P1.muE/(P1.RE+Y(3))^2;
dY    = [(F_S-F_D)/Y(5) - g*sin(Y(2));
        -(g/Y(1)-Y(1)/(P1.RE+Y(3)))*cos(Y(2));
         Y(1)*sin(Y(2));
         Y(1)*cos(Y(2));
         -0];
end

function [value,isterminal,direction] = myeventGamma(t, Y)
        value = (Y(2)) ; % detect gamma = 0
        isterminal = 1;  % stop the integration
        direction = 0;   % negative direction
end


function [value,isterminal,direction] = myeventNick(t, Y)
        value = (Y(3)-900) ;   % detect h = 900 m
        isterminal = 1;  % stop the integration
        direction = 0;   % negative direction
end

function [value,isterminal,direction] = myeventxend(t, Y)
        value = (Y(4)-1.2e6) ;   % detect x = 1200 km
        isterminal = 1;  % stop the integration
        direction = 0;   % negative direction
end
