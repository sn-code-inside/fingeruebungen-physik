% -------------------------------------------------------------------------
% ReentryReal.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Dynamik des Wiedereintritts eines Raumschiffes
% in die Erdatmosphäre.
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
RE  = 6.378e3;                     % Erdradiuse in km
muE = G*ME*1e-9;                   % µ in km3 /kg /s2
gE  = muE/RE^2;                    % gE in km/s^2
hA  = 400;                         % Anfangshöhe in km
hRe = 122;                         % Reentry Höhe in km
rA      = RE + hA;
rRe     = RE + hRe;
vA      = sqrt(muE/rA);
q       = rA/rRe;
qminus  = q-1;


%% Anfangswerte 

gammaRe = deg2rad(5);                             % FPA @ Reentry
vRe     = vA*(1+0.75*qminus-gammaRe^2/8/qminus);  % v @ Reentry

Deltav  = 0.25*vA*(qminus+gammaRe^2/2/qminus);    % Delta v for Reentry
% oder genauer
Deltav  = cos(gammaRe)*vRe/q - vA;                % Delta v for Reentry
exz     = (q^2-(2*q-1)*(cos(gammaRe))^2)/(q^2-(cos(gammaRe))^2);
a       = rA/(1+exz);
rP      = a*(1-exz);
vARe    = sqrt(muE*(2/rA-1/a));
upsRe   = 180-acosd((qminus^2*(cos(gammaRe))^2-q^2*(sin(gammaRe))^2)/...
                (qminus^2+q^2*(sin(gammaRe))^2));
v0      = vRe;
gamma0  = gammaRe;
h0      = hRe;
x0      = 0;

% Parameter der DGL
mL         = 5e03;                % Payload in kg
rho0       = 1.752e09;            % Dichte in kg/km^3
H0         = 6.7;                 % H0 in km
kappaD     = 25;                  % normierte Widerstandswert
cwD        = 2*kappaD*mL/rho0/H0; % Drag Widerstandsbeiwert x Fläche km^2 
facAD      = [0;0.1;0.2;0.3;0.4;0.5];
cwAvec     = facAD*cwD;           % Auftriebswiderstandswert x Fläche km^2 

% Barometrische Höhenformel
rho        = @(h)rho0*exp(-h/H0);  

% AB fÃ¼r DGL
P1.mL     = mL;
P1.gE     = gE;
P1.RE     = RE;
P1.cwD    = cwD;

%% MATLABs Runge-Kutta ode45 Routine 
opts   = odeset('AbsTol',1.e-9,'RelTol',1.e-7);
AB1    = [v0;gamma0;h0;x0];          % AB fÃ¼r DGL
tmax = round(1.25*sqrt(2*hRe/gE/sin(gammaRe))/100)*100;
xmax = round(1.5*hRe/tan(gammaRe)/200)*200;
hmax = 125;
tspan1 = linspace(0,tmax,1001);
% Schleife über verschiede Drag/Auftriebsverhältnisse
for k=1:length(cwAvec)
    P1.facAD    = facAD(k);
    [t,Y1]      = ode45(@DGL_Reentry, tspan1, AB1, opts, P1, rho);
    v1(k,:)     = Y1(:,1);
    gamma1(k,:) = Y1(:,2);
    x1(k,:)     = Y1(:,4);
    h1(k,:)     = Y1(:,3);
    t1(k,:)     = t(:);
end  
% Berechnung Beschleunigungen
for k=1:length(cwAvec)
    a1(k,:) = zeros(1,length(v1(k,:)));
    for m=2:length(v1(k,:))
        a1(k,m) = (v1(k,m)-v1(k,m-1))/(t1(k,m)-t1(k,m-1))/gE;
    end
    a1(k,1) = a1(k,2);
end

%% Graphische Ausgabe

strgamma = string(num2str(rad2deg(gammaRe),'%3.0f'));
figure('name','h-x-Plot')
hold on
for k=1:length(cwAvec)
    p(k)=plot(x1(k,:),h1(k,:),'color',Colors(k,:),'Linewidth',2);
    lgdstr(k,:) = strjoin(["A/D = ",...
                  string(num2str(facAD(k),'%3.1f'))]);
end    
grid on
str= strjoin(["h(x),  \gamma_{Re} := ",strgamma,'°']);
grid on;
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'Location','northeast','NumColumns',2); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, xmax, 0 hmax])
xlabel('x in km','FontSize',14); 
ylabel('h in km','FontSize',14)
legend box off;
set(gca,'FontSize',16);

strgamma = string(num2str(rad2deg(gammaRe),'%3.0f'));
figure('name','v-h-Plot')
hold on
for k=1:length(cwAvec)
    p(k)=plot(h1(k,:),v1(k,:),'color',Colors(k,:),'Linewidth',2);
    lgdstr(k,:) = strjoin(["A/D = ",...
                  string(num2str(facAD(k),'%3.1f'))]);
end    
grid on
str= strjoin(["v(h),  \gamma_{Re} := ",strgamma,'°']);
grid on;
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'Location','northeast','NumColumns',2); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, hmax, 0 1.1*vRe]);
xlabel('h in km','FontSize',14); 
ylabel('v in km/s','FontSize',14)
legend box off;
set(gca, 'XDir','reverse')
set(gca,'FontSize',16);

figure('name','h-t-Plot')
subplot(2,1,1)
maxv  = 1.2*max(max(v1));
hold on
for k=1:length(cwAvec)
    p(k)= plot(t1(k,:),h1(k,:),'color',Colors(k,:),'Linewidth',2);
end    
str= strjoin(["h(t), \gamma_{Re} := ",strgamma,'°']);
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'Location','northeast','NumColumns',2); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, tmax, 0 hmax]);
legend box off;
xlabel('t in s','FontSize',14); 
ylabel('h in km','FontSize',14)
grid on
set(gca,'FontSize',16);

dde23

subplot(2,1,2)
hold on
for k=1:length(cwAvec)
    p(k)= plot(t1(k,:),v1(k,:),'color',Colors(k,:),'Linewidth',2);
end    
str= strjoin(["v(t), \gamma_{Re} := ",strgamma,'°']);
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'Location','northeast','NumColumns',2); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, tmax, 0 round(1.5*vRe)]);
xlabel('t in s','FontSize',14); 
ylabel('v in km/s','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);

%---------------------------------------------------------------------
figure('name','a(t)-gamma(t)-Plot')
subplot(2,1,1)
maxa  = 1.2*max(max(a1));
mina  = 1.2*min(min(a1));
hold on
for k=1:length(cwAvec)
    p(k)= plot(t1(k,:),a1(k,:),'color',Colors(k,:),'Linewidth',2);
end    
str= strjoin(["a(t), \gamma_{Re} := ",strgamma,'°']);
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'Location','southeast','NumColumns',2); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, tmax, mina maxa]);
xlabel('t in s','FontSize',14); 
ylabel('Beschleunigung x g','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);

subplot(2,1,2)
maxangle  = max(gamma1(length(cwAvec),:));
minangle  = min(min(gamma1(:,:)));
maxwinkel = 10*round(rad2deg(maxangle)/10);
minwinkel = 20*round(rad2deg(minangle)/20)-10;
hold on
for k=1:length(cwAvec)
    p(k)= plot(t1(k,:),rad2deg(gamma1(k,:)),...
          'color',Colors(k,:),'Linewidth',2);
end    
str= strjoin(["FPA(t), \gamma_{Re} := ",strgamma,'°']);
hp1 = title(str,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:k),lgdstr(1:k,:),'location','southeast','NumColumns',2);
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, tmax, -20 100]);
xlabel('t in s','FontSize',14); 
ylabel('\gamma in °','FontSize',14)
legend box off;
grid on
set(gca,'FontSize',16);


%%

% DGL
function dY = DGL_Reentry(t, Y, P1, Dichte)
% Y(1):v, Y(2):gamma, Y(3):h, Y(4):x 
mR    = P1.mL;
F_D   = 0.5*feval(Dichte,Y(3))*P1.cwD*Y(1)^2;
F_A   = F_D*P1.facAD;
g     = P1.gE./(1+2*Y(3)./P1.RE+(Y(3)/P1.RE)^2);
dY    = [-F_D/mR + g*sin(Y(2));
         (-F_A/mR/Y(1) +(g/Y(1)-Y(1)/(P1.RE+Y(3)))*cos(Y(2)));
         - Y(1)*sin(Y(2));
         Y(1)*cos(Y(2))];
end

