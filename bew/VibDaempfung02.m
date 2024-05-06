% -------------------------------------------------------------------------
% VibDaempfung02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Schwingungsdämpfung einer Deckenaufhängung
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Numerik
% Parameter

M       =   50;            % Masse Gerät in kg
m       =   25;            % Masse Gegenmasse in kg
kM      =  15000;          % Federkonstante Gerät
kT      =  25000;          % Federkonstante Gegenmasse
etaM    =  150;            % Dämpfungskoeffizient
etaT    =  30;            % Dämpfungskoeffizient
tmax    =   3;

% M       =   1;            % Masse Gerät in kg
% m       =  0.6;           % Masse Gegenmasse in kg
% kM      =  120;           % Federkonstante Gerät
% kT      =  135;           % Federkonstante Gegenmasse
% etaM    =  1.4;           % Dämpfungskoeffizient
% etaT    =  1.0;           % Dämpfungskoeffizient
% tmax    =  5;

Ax    = 0.10;            % Anfangsamplitude externe Anregung;
omegaM  = sqrt(kM/M);
omegaT  = sqrt(kT/m);
Omega   = omegaT;

fprintf('\n ');
fprintf('\n ');
fprintf('\n Masse Mikroskop :            M    = %8.2f kg', M);
fprintf('\n Masse Gegenmasse :           m    = %8.2f kg', m);
fprintf('\n Federkonstante Gerät :       kM   = %8.2f kg/s² (N/m)', kM);
fprintf('\n Federkonstante Gegenmasse :  kT   = %8.2f kg/s² (N/m)', kT);
fprintf('\n Dämpfung Gerät :             etaM = %8.2f kg/s  (N/(m/s)', etaM);
fprintf('\n Dämpfung Gegenmasse :        etaT = %8.2f kg/s  (N/(m/s)', etaT);
fprintf('\n Eigenfrequenz Gerät :        fM   = %8.2f 1/s  ', omegaM/2/pi);
fprintf('\n Eigenfrequenz Gegenmasse :   fT   = %8.2f 1/s  ', omegaT/2/pi);
fprintf('\n Anregungsfrequenz  :         fext = %8.2f 1/s  ', Omega/2/pi);
fprintf('\n Anregungsamplitude  :        Ax   = %8.2f 1/s  ', Ax);
fprintf('\n ');

% Parameter in LGL
% Anfangswerte
zM0  = Ax;
dzM0 = 0;
zT0  = 0;
dzT0 = 0;

%% Berechnungen ODE45

%Anfangswerte
AB=[zM0;dzM0;zT0;dzT0]; % AB für ode45

P1.M     = M;
P1.m     = m;
P1.kT    = kT;
P1.kM    = kM;
P1.etaT  = etaT;
P1.etaM  = etaM;
P1.Ax    = Ax;
P1.Omega = Omega;
omega = linspace(2*pi*0.1,2*pi*10,500);

gamma1 = atan2(etaT*omega,(kT-m*omega.^2));
gamma2 = atan2(etaM*omega,(kT+kM-M*omega.^2));
Q1     = sqrt((kM+kT-M*omega.*omega).^2+(etaM*omega).^2);
Q2     = kT^2./ sqrt((kT-m*omega.*omega).^2+(etaT*omega).^2);


% Amplituden
AM     = kM*Ax./sqrt(Q1.^2+Q2.^2-2*Q1.*Q2.*cos(gamma1+gamma2));
AT     = kT.*AM./sqrt((kT-m*omega.^2).^2+(etaT*omega).^2);
% Ohne Tilgung
A0     = kM*Ax/M./sqrt((kM/M-omega.^2).^2+etaM^2/M^2.*omega.^2);

% Numerische Lösung Lagrangegleichung


opt=odeset('AbsTol',1.e-7,'RelTol',1.e-7);
% Numerische Lösung LGL mit Reibung
[t, Y] = ode45(@dgl_VibDamp02,[0.0,tmax],AB,opt,P1); 
zM  = Y(:,1);
dzM = Y(:,2); 
zT  = Y(:,3);

P1.kT =0;  % keine Kopplung
[t0, Y0] = ode45(@dgl_VibDamp02,[0.0,tmax],AB,opt,P1); 
zMung  = Y0(:,1);


f = omega/2/pi;
%% 
% Graphische Ausgabe

% Frequenzgang 
figure();
p(1) = semilogy(f, AM,'Color',Colors(2,:), 'LineWidth',2);
hold on
p(2) = semilogy(f, AT,'Color',Colors(3,:), 'LineWidth',1,...
       'LineStyle', Style(2));
p(3) = semilogy(f, A0,'Color',Colors(4,:), 'LineWidth',1,...
       'LineStyle', Style(3));
grid on
ylabel('Amplitude in m','FontSize',14)
axis ([0,omegaT/pi,0.001,1]);
xlabel('Frequenz f in Hz','FontSize',14)
h=title('Deckenaufhängung ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p,'mit Tilgung ','Gegengewicht ', 'ohne Tilgung','location','southwest',...
         'NumColumns',1);
legend box off
set(gca,'FontSize',16);


% Dynamik 1
figure();
p(1) = plot(t, zM,'Color',Colors(2,:), 'LineWidth',2);
hold on
p(2) = plot(t, zT,'Color',Colors(3,:), 'LineWidth',1,...
    'LineStyle', Style(2));
p(3) = plot(t0, zMung,'Color',Colors(4,:), 'LineWidth',1,...
    'LineStyle', Style(3));
axis ([0,tmax,-inf,inf]);
grid on
ylabel('Amplitude in m','FontSize',14)
xlabel('Zeit in s','FontSize',14)
h=title('Deckenaufhängung ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p,'mit Tilgung ','Gegengewicht ', 'ohne Tilgung','location','southeast',...
         'NumColumns',1);
legend box off
set(gca,'FontSize',16);

% Dynamik 2
figure();
p(1) = semilogy(t, abs(zM),'Color',Colors(2,:), 'LineWidth',2);
hold on
p(2) = semilogy(t, abs(zT),'Color',Colors(3,:), 'LineWidth',1,...
    'LineStyle', Style(2));
p(2) = semilogy(t0, abs(zMung),'Color',Colors(4,:), 'LineWidth',1,...
    'LineStyle', Style(3));
axis ([0,tmax,0.001,1]);
grid on
ylabel('Amplitude in m','FontSize',14)
xlabel('Zeit in s','FontSize',14)
h=title('Deckenaufhängung ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p,'mit Tilgung ','Gegengewicht ', 'ohne Tilgung','location','northeast',...
         'NumColumns',1);
legend box off
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%% Funktionen

% Lagrangegleichung
function dY = dgl_VibDamp02(t,Y,P1)
    % Y(1)- Position zM(t), 
    % Y(2)- Geschwindigkeit dotzM(t)
    % Y(3)- Position zm(t), 
    % Y(4)- Geschwindigkeit dotzM(t)
    zM  = Y(1);
    dzM = Y(2);
    zT  = Y(3);
    dzT = Y(4);
    M    = P1.M;
    m    = P1.m;
    kM   = P1.kM;
    kT   = P1.kT;
    etaM = P1.etaM; 
    etaT = P1.etaT; 
    Ax   = P1.Ax;
    Omega= P1.Omega;
    dY = [dzM;...
      -(etaM*dzM + (kM+kT)*zM - kT*zT + kM*Ax*cos(Omega*t))/M;
      dzT;
      -(etaT*dzT+kT*(zT-zM))/m];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

