% -------------------------------------------------------------------------
% RaumschiffStartReal.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet Start einer einstufigen Rakete von der Erdoberfläche
% unter Berücksichtigung von Reibung und Nickwinkel-Korrekturen.
%
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
G   = 6.67430e-11;                 % G in m3 /kg /s2
ME  = 5.98e24;                     % Masse der Erde in kg
RE  = 6.378e6;                     % Erdradiuse in m
gE  = G*ME/RE^2*1e-3;              % gE in km/s^2

% Anfangswerte 
v0     = 0.0001;
gamma0 = pi/2;
h0     = 0;
x0     = 0;

% Parameter der DGL

Del_gamma  = deg2rad([0,1,2.5,5,7.5]);          % Nickkorrektur
cwD        = 1.5e-6;              % Widerstandsbeiwert x Fläche in km^2 
m0         = 1e5;                 % Startmasse in kg
mL         = 5e3;                 % Payload in kg
dm         = 500;                 % Masseverlust in kg/s
ISpez      = 300;                 % spez. Impuls in s  Ispez = FS/g/dm;
FS         = ISpez*gE*dm;         % Schub in N
hNick      = 1000;                % Höhe der Nickkorrektur
tBurnOut   = (m0-mL)/dm;          % Burn-out in s
rho0       = 1.752e09;            % Dichte in kg/km^3
H0         = 6.7;                 % H0 in km

% Barometrische Höhenformel
rho        = @(h)rho0*exp(-h/H0);  

AB = [v0;gamma0;h0;x0];          % AB für DGL
P1.gE     = gE;
P1.RE     = RE;
P1.cwD    = cwD;
P1.m0     = m0;
P1.mf     = mL;
P1.dm     = dm;
P1.tBurnOut = tBurnOut;
P1.FS       = FS;

% MATLABs Runge-Kutta ode45 Routine 

opts   = odeset('AbsTol',1.e-9,'RelTol',1.e-7);
AB1    = [v0;gamma0;h0;x0];          % AB für DGL
tspan0 = linspace(0,tBurnOut+40,251);
% Ohne Reibung 
[t0,Y0]= ode45(@DGL_Launch_ZeroDrag,tspan0, AB1,opts,P1,rho);
%Anfangsphase (t < 20s bis Nickkorrektur
tspan1 = linspace(0,20,251);
[t1,Y1]= ode45(@DGL_Launch,tspan1, AB1,opts,P1,rho);
tend1  = t1(end);
%Phase bis Burn-Out, nach Nickkorrektur
tspan2 = linspace(tend1,tBurnOut,200);
for k=1:length(Del_gamma)
    AB2    = [Y1(end,1);Y1(end,2)-Del_gamma(k);Y1(end,3);Y1(end,4)]; 
    [t,Y2]= ode45(@DGL_Launch, tspan2, AB2, opts, P1, rho);
    v2(k,:)     = Y2(:,1);
    gamma2(k,:) = Y2(:,2);
    x2(k,:)     = Y2(:,4);
    h2(k,:)     = Y2(:,3);
    t2(k,:)     = t(:);
end    
tend2  = tBurnOut;
tspan3 = linspace(tend2,tend2+40,101);
for k=1:length(Del_gamma)
    AB2    = [v2(k,end); gamma2(k,end);h2(k,end);x2(k,end)]; 
    [t,Y3]= ode45(@DGL_Launch, tspan3, AB2, opts, P1, rho);
    v3(k,:)     = Y3(:,1);
    gamma3(k,:) = Y3(:,2);
    x3(k,:)     = Y3(:,4);
    h3(k,:)     = Y3(:,3);
    t3(k,:)     = t(:);
end    


%% Graphische Ausgabe

figure('name','h-x-Plot')
hold on
plot(Y1(:,4),Y1(:,3),'color',Colors(1,:));
for k=1:length(Del_gamma)
    p(k)=plot(x2(k,:),h2(k,:),'color',Colors(k,:),'Linewidth',2);
    plot(x3(k,:),h3(k,:),'color',Colors(k,:),...
           'Linewidth',2,'LineStyle',':');
    lgdstr(k,:) = strjoin(["\Delta \gamma = ",...
                  string(num2str(rad2deg(Del_gamma(k)),'%3.1f'))]);
end    
p(k+1)=plot(Y0(:,4)-25,Y0(:,3),'color',Colors(8,:),'Linewidth',1,...
              'LineStyle','-');
lgdstr(k+1,:) = "ohne Drag";
grid on
str= "h(x)-Graphik";
grid on;
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([-100, 600, 0 600]);
xlabel('x in km','FontSize',14); 
ylabel('h in km','FontSize',14)
legend box off;
set(gca,'FontSize',16);



figure('name','h-t-Plot')
hold on
plot(t1(:),Y1(:,3),'color',Colors(1,:));
for k=1:length(Del_gamma)
    p(k)= plot(t2(k,:),h2(k,:),'color',Colors(k,:),'Linewidth',2);
    plot(t3(k,:),h3(k,:),'color',Colors(k,:),...
           'Linewidth',2,'LineStyle',':');
    lgdstr(k,:) = strjoin(["\Delta \gamma = ",...
                  string(num2str(rad2deg(Del_gamma(k)),'%3.1f'))]);
end    
p(k+1)=plot(t0(:),Y0(:,3),'color',Colors(8,:),'Linewidth',1,...
              'LineStyle','-');
str= "h(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'location','northwest'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, 300, 0 600]);
legend box off;
xlabel('t in s','FontSize',14); 
ylabel('h in km','FontSize',14)
grid on
set(gca,'FontSize',16);

figure('name','v(t)-Plot')
subplot(2,1,1)
hold on
plot(t1(:),Y1(:,1),'color',Colors(1,:));
for k=1:length(Del_gamma)
    p(k)= plot(t2(k,:),v2(k,:),'color',Colors(k,:),'Linewidth',2);
    plot(t3(k,:),v3(k,:),'color',Colors(k,:),...
           'Linewidth',2,'LineStyle',':');
end    
p(k+1)=plot(t0(:),Y0(:,1),'color',Colors(8,:),'Linewidth',1,...
              'LineStyle','-');
str= "v(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'location','northwest'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, 300, 0 10]);
xlabel('t in s','FontSize',14); 
ylabel('v in km/s','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);

subplot(2,1,2)
hold on
plot(t1(:),rad2deg(Y1(:,2)),'color',Colors(1,:));
for k=1:length(Del_gamma)
    p(k)= plot(t2(k,:),rad2deg(gamma2(k,:)),'color',Colors(k,:),'Linewidth',2);
    plot(t3(k,:),rad2deg(gamma3(k,:)),'color',Colors(k,:),...
           'Linewidth',2,'LineStyle',':');
end    
str= "\gamma(t) - Graphik";
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:k),lgdstr(1:k,:),'location','southwest'); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, 300, 30 90]);
xlabel('t in s','FontSize',14); 
ylabel('\gamma in °','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);







%%

% DGL
function dY = DGL_Launch(t, Y, P1, Dichte)
% Y(1):v, Y(2):gamma, Y(3):h, Y(4):x 
F_S   = P1.FS*(t<P1.tBurnOut);  
mR    = (P1.m0 - P1.dm*t)*(t<P1.tBurnOut)+P1.mf*(t>=P1.tBurnOut);
F_D   = 0.5*feval(Dichte,Y(3))*P1.cwD*Y(1)^2;
g     = P1.gE/(1+2*Y(3)/P1.RE+(Y(3)/P1.RE)^2);
dY    = [(F_S-F_D)/mR - g*sin(Y(2));
        -(g/Y(1)-Y(1)/(P1.RE+Y(3)))*cos(Y(2));
         Y(1)*sin(Y(2));
         Y(1)*cos(Y(2))];
end

function dY = DGL_Launch_ZeroDrag(t, Y, P1, Dichte)
% Y(1):v, Y(2):gamma, Y(3):h, Y(4):x 
F_S   = P1.FS*(t<P1.tBurnOut);  
mR    = (P1.m0 - P1.dm*t)*(t<P1.tBurnOut)+P1.mf*(t>=P1.tBurnOut);
F_D   = 0;
g     = P1.gE/(1+2*Y(3)/P1.RE+(Y(3)/P1.RE)^2);
dY    = [(F_S-F_D)/mR - g*sin(Y(2));
        -(g/Y(1)-Y(1)/(P1.RE+Y(3)))*cos(Y(2));
         Y(1)*sin(Y(2));
         Y(1)*cos(Y(2))];
end
