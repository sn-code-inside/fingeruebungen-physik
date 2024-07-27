% -------------------------------------------------------------------------
% RakentenStart2Orbit02.m
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
% Variabel: Raketenparameter, Nutzlast, Resttreibstoff
% Bsp. Ariane 1
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
Del_gamma_v = deg2rad([0.5,0.8,1.5,2.5]);  % mögliche Nickkorrekturen
Del_gamma   = Del_gamma_v(2);              % gewählte Nickkorrektur
fracBurn    = 0.75;                        % Anteil Treibstoffuse Stufe 3

Parastr(1,:)= sprintf('Nutzlast : %6.1f kg', mLA);
Parastr(2,:)= sprintf('Nickkorr.: %6.1f °', rad2deg(Del_gamma));
Parastr(3,:)= sprintf('Used Fuel: %6.3f   ', fracBurn);

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

% Anfangswerte 
v0     = 0.001;
gamma0 = pi/2;
h0     = 0;
x0     = 0;
m0     = m0A(1);

% MATLABs Runge-Kutta ode45 Routine 
opts   = odeset('AbsTol',1.e-09,'RelTol',1.e-08);

% Berechnung Anfangsphase (bis Nickkorrektur)
AB0       = [v0;gamma0;h0;x0;m0];          % AB fÃ¼r DGL
tNick     = 20;                            % Zeit der Nickkorrektur in s
P1.Stage  = 1;
tspan1    = linspace(0,tNick,101);
[t1,Y]    = ode45(@DGL_Thrust,tspan1, AB0,opts,P1);
v1(:)     = Y(:,1);
gamma1(:) = Y(:,2);
h1(:)     = Y(:,3);
x1(:)     = Y(:,4);
m1(:)     = Y(:,5);

% Berechnung erste Stufe bis Burn-Out Stufe 1

AB1       = [v1(end);gamma1(end)-Del_gamma;h1(end);x1(end);m1(end)];    
P1.Stage  = 1;
tspan1    = linspace(tNick,tBA(1),501);
[t2,Y]    = ode45(@DGL_Thrust,tspan1, AB1,opts,P1);
v2(:)     = Y(:,1);
gamma2(:) = Y(:,2);
h2(:)     = Y(:,3);
x2(:)     = Y(:,4);
m2(:)     = Y(:,5);

% Berechnung zweite Stufe bis Burn-Out Stufe 2
AB2       = [v2(end);gamma2(end);h2(end);x2(end);m0A(2)];    
P1.Stage  = 2;
tspan2    = linspace(t2(end),t2(end)+tBA(2),501);
[t3,Y]    = ode45(@DGL_Thrust,tspan2, AB2,opts,P1);
v3(:)     = Y(:,1);
gamma3(:) = Y(:,2);
h3(:)     = Y(:,3);
x3(:)     = Y(:,4);
m3(:)     = Y(:,5);

% Berechnung dritte Stufe bis Burn-Out Stufe 3
AB3       = [v3(end);gamma3(end);h3(end);x3(end);m0A(3)];    
% opts3     = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'Events',@myevent);
P1.Stage  = 3;
tspan3    = linspace(t3(end),t3(end)+tBA(3)*fracBurn,501);
% [t3,Y,te,ye,ie] = ode45(@(t3,Y) DGL_Launch(t3,Y,P1), tspan3, AB3, opts3);
[t4,Y]    = ode45(@DGL_Thrust, tspan3, AB3,opts,P1);
v4(:)     = Y(:,1);
gamma4(:) = Y(:,2);
h4(:)     = Y(:,3);
x4(:)     = Y(:,4);
m4(:)     = Y(:,5);

 
%Phase nach Burn-Out, free trajectory
if fracBurn == 1
    mfinal = mLA;
else
    mfinal = m4(end);
end
AB4       = [v4(end);gamma4(end);h4(end);x4(end);mfinal];   
tspan4    = linspace(t4(end),t4(end)+300,501);
[t5,Y]    = ode45(@DGL_FreeTraj, tspan4, AB4,opts,P1);
v5(:)     = Y(:,1);
gamma5(:) = Y(:,2);
x5(:)     = Y(:,4);
h5(:)     = Y(:,3);
m5(:)     = Y(:,5);

%% Graphische Ausgabe
%--------------------------------------------------------------------------

figure('name','h-x-Plot Ariane-1')
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
axis([0, inf, 0 1.2*inf]);
xlabel('x in km','FontSize',14); 
ylabel('h in km','FontSize',14)
for k=1:3
  PP=text(1000,100-k*25,Parastr(k,:));
  set(PP,'FontSize',14);
end
legend box off;
set(gca,'FontSize',16);

figure('name','Übersicht-Plots Ariane-1')
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
axis([0, t5(end), 0 1.2*max(h4)/1000]);
for k=1:3
  PP=text(25,max(h4)/1000-k*35,Parastr(k,:));
  set(PP,'FontSize',14);
end
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
axis([0, t5(end), 0 1.2*max(v4/1000)]);
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
ylabel('\gamma in °','FontSize',14)
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
hp2=legend(p(1:k),lgdstr(1:k,:),'location','northeast'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, t5(end), 1e3 2*m0]);
xlabel('t in s','FontSize',14); 
ylabel('m in kg','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);


%%  Funktionen
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
F_D   = 0.5*rho*P1.cwD*Y(1)^2;
g     = P1.muE/(P1.RE+Y(3))^2;
dY    = [-F_D/Y(5) - g*sin(Y(2));
        -(g/Y(1)-Y(1)/(P1.RE+Y(3)))*cos(Y(2));
         Y(1)*sin(Y(2));
         Y(1)*cos(Y(2));
         -0];
end

function [value,isterminal,direction] = myevent(t, Y)
        value = (Y(2)) ; % detect gamma = 0
        isterminal = 1;  % stop the integration
        direction = 0;   % negative direction
end

function [value,isterminal,direction] = myeventxend(t, Y)
        value = (Y(4)-3e6) ;   % detect x = 1200 km
        isterminal = 1;  % stop the integration
        direction = 0;   % negative direction
end
