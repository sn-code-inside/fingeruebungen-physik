% -------------------------------------------------------------------------
% DispersionWelle.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm zeigt Dispersion am Unterschied Gruppen- und Phasen-
% geschwindigkeit.
% 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
%normal Dispersion  (Auswahl)
lambda1  = 0.95;                        % Wellenlänge 1 in m
lambda2  = 1.05;                        % Wellenlänge 2 in m
lambda   = (lambda1+lambda2)/2;         % mittlere Wellenlänge in m
vph1     = 0.975;                       % Phasengeschwindigkeit 1 in m/s
vph2     = 1.025;                       % Phasengeschwindigkeit 1 in m/s

%anormale Dispersion (Auswahl)
% lambda1  = 0.95;                        % Wellenlänge 1 in m
% lambda2  = 1.05;                        % Wellenlänge 2 in m
% lambda   = (lambda1+lambda2)/2;         % mittlere Wellenlänge in m
% vph1     = 1.025;                       % Phasengeschwindigkeit 1 in m/s
% vph2     = 0.975;                       % Phasengeschwindigkeit 1 in m/s

k0       = 2*pi/lambda ;                % Wellenzahl  in m/s
k1       = 2*pi/lambda1;                % Wellenzahl 1 in m/s
k2       = 2*pi/lambda2;                % Wellenzahl 2 in m/s
omega1   = vph1*k1;                       % Kreisfrequenz 1 in 1/s
omega2   = vph2*k2;                       % Kreisfrequenz 2 in 1/s

vph      = (vph1+vph2)/2;               % Näherung mittl. Phasengeschwindigkeit
vgr      = (omega1-omega2)/(k1-k2); 
% Näherung Gruppengeschwindigkeit
vgr2     = vph-lambda*(vph1-vph2)/(lambda1-lambda2); 
% Näherung GruppenPhasengeschwindigkeit

xend = 20;
xsteps =400;
x = linspace(0,xend,xsteps);
deltax = xend/xsteps;

tend   = 40;
tsteps = 40;
t = linspace(0,tend,tsteps);
deltat = tend/tsteps;

xmax  = 20;
ymax  = 2.25;


%% 1. Zwei Wellen in gleiche Richtung laufend

% Berechnung Wellen
% Einzelwellen
psi1 = zeros(length(x),length(t));
psi2 = zeros(length(x),length(t));
% Mittlere Ausbreitung (Trägerwelle)
psi  = zeros(length(x),length(t));
for k=1:length(t)
  for m=1:length(x)
    psi (m,k) = cos(k0*x(m)/lambda-vph*k0*t(k));
    psi1(m,k) = cos(k1*x(m)-omega1*t(k));
    psi2(m,k) = cos(k2*x(m)-omega2*t(k));
    psig(m,k) = psi1(m,k)+psi2(m,k);
  end
end


% Simulation zeitlicher Verrlauf
% figure('name','Simulation 2 Wellen')
% for k=1:length(t)
%      clf
%      hold on 
%      plot(x(:), psi(:,k));
%      plot(x(:), psig(:,k));
%      grid on;
%      pause(0.01);
%      hold off
% end
% grid on;
% axis([0 xmax -ymax ymax]);
% xlabel('{\itx}');
% ylabel('{\it\psi(\itx, \itt})');
% set(gca,'FontSize',16,'FontName','Times');


% Wellenbild zu zwei verschiedenen Zeiten
figure('name','Wellenbild zu zwei verschiedenen Zeiten')
subplot(2,1,1)
hold on 
tmess = 0;
kmess = 1;
plot(x(:), psi1(:,kmess),'color',Colors(2,:),'linewidth',1);
plot(x(:), psig(:,kmess),'color',Colors(4,:));
[up,lo] = envelope(psig(:,kmess),xsteps,'analytic');
plot(x,up,'color',Colors(4,:),'linewidth',1);
plot(x,lo,'color',Colors(4,:),'linewidth',1);
line([0.01 0.01],[-2 2],'color',Colors(4,:),'linewidth',2);
line([0.05 0.05],[-1 1],'color',Colors(2,:),'linewidth',2);
grid on;
tmess_str = strcat('t = ',num2str(tmess,3));
tmess_str = strcat(tmess_str,' s');
axis([0 xmax -ymax ymax]);
xlabel('{\itx}');
ylabel('{\it\psi(\itx, \itt})');
h2 = title(tmess_str);
set(h2, 'FontSize',14,'FontName','Times', 'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');

tmess = 10;
kmess = tmess/deltat;
subplot(2,1,2)
hold on 
plot(x(:), psi1(:,kmess),'color',Colors(2,:));
plot(x(:), psig(:,kmess),'color',Colors(4,:));
[up,lo] = envelope(psig(:,kmess),xsteps,'analytic');
plot(x,up,'color',Colors(4,:));
plot(x,lo,'color',Colors(4,:));
[~,mup]=max(up);
line([mup*deltax mup*deltax],[-2 2],'color',Colors(4,:),'linewidth',2);
line([vph/vgr*mup*deltax vph/vgr*mup*deltax],[-1 1],'color',...
     Colors(4,:),'linewidth',1);
line([vph1/vgr*mup*deltax vph1/vgr*mup*deltax],[-1 1],'color',...
     Colors(2,:),'linewidth',2);
grid on;
tmess_str = strcat('t = ',num2str(tmess,3));
tmess_str = strcat(tmess_str,' s');
axis([0 xmax -ymax ymax]);
xlabel('{\itx}');
ylabel('{\it\psi(\itx})');
h2 = title(tmess_str);
set(h2, 'FontSize',14,'FontName','Times', 'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');



%% 3. Wellenpaket

% Berechnung Wellenpaket aus Einzelwellen
psi_wp  = zeros(length(x),length(t));

vph0    = 1;
vd      = 0;                            % Dispersionsparameter = 0 m/s
vd      = -0.4;                         % Dispersionsparameter = -0.7 m/s
lambda0 = 1.0;                          % Wellenlänge 1 in m
kw0     = 2*pi/lambda0 ;                % Wellenzahl  in m/s
omega0  = vph0*kw0;                     % Kreisfrequenz 1 in 1/s
Nmax    = 101;
% Mittlere Ausbreitung (Trägerwelle)
for n=1:Nmax
    lambda(n) = lambda0 + 0.2*(n-Nmax)/(Nmax-1); % Wellenlänge
    kw(n)     = 2*pi/lambda(n);
    vph(n)    = vph0 + vd*(abs(kw(n))-kw0)/kw0;
    omega(n)  =  vph(n)*kw(n);
    A(n)   = exp(-((n-Nmax)/(Nmax-1))^2/4);
    for k=1:length(t)
        for m=1:length(x)
            psi_wp(m,k) = psi_wp(m,k) + A(n)*cos(kw(n)*x(m)-omega(n)*t(k));
        end
    end
end

%Simulation zeitlicher Verrlauf
figure('name','Simulation Ausbreitung Wellenpaket')
xmax = 20;
ymax = 100;
axis([0 xmax -ymax ymax]);
for k=1:length(t)
    clf
    hold on
    plot(x,psi_wp(:,k));
    axis([0 xmax -ymax ymax]);
    grid on;
    pause(0.04);
    hold off
end
grid on;
axis([0 xmax -ymax ymax]);
xlabel('{\itx}');
ylabel('{\it\psi(\itx, \itt})');
set(gca,'FontSize',16,'FontName','Times');


% Wellenbild zu zwei verschiedenen Zeiten
figure('name','Wellenpaket zu zwei verschiedenen Zeiten')
subplot(2,1,1)
hold on 
tmess = 0;
kmess = 1;
plot(x(:), psi_wp(:,kmess),'color',Colors(4,:));
[up,lo] = envelope(psi_wp(:,kmess),xsteps,'analytic');
plot(x,up,'color',Colors(4,:),'linewidth',1);
plot(x,lo,'color',Colors(4,:),'linewidth',1);
grid on;
tmess_str = strcat('t = ',num2str(tmess,3));
tmess_str = strcat(tmess_str,' s');
axis([0 xmax -ymax ymax]);
xlabel('{\itx}');
ylabel('{\it\psi(\itx, \itt})');
h2 = title(tmess_str);
set(h2, 'FontSize',14,'FontName','Times', 'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');

tmess = 30;
kmess = tmess/deltat;
subplot(2,1,2)
hold on 
plot(x(:), psi_wp(:,kmess),'color',Colors(4,:));
[up,lo] = envelope(psi_wp(:,kmess),xsteps,'analytic');
plot(x,up,'color',Colors(4,:),'linewidth',1);
plot(x,lo,'color',Colors(4,:),'linewidth',1);
plot(x,lo,'color',Colors(4,:));
tmess_str = strcat('t = ',num2str(tmess,3));
tmess_str = strcat(tmess_str,' s');
grid on
axis([0 xmax -ymax ymax]);
xlabel('{\itx}');
ylabel('{\it\psi(\itx})');
h2 = title(tmess_str);
set(h2, 'FontSize',14,'FontName','Times', 'FontWeight','normal');
set(gca,'FontSize',14,'FontName','Times');


% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

