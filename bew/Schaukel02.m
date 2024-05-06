% -------------------------------------------------------------------------
% Schaukel02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Schaukel als Getriebener Oszillator
% 
% Lösung der DGL für die Schaukel als getriebener Oszillator
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Definiert die verallgemeinerten Koordinaten theta1, theta2 und Parameter

syms J1 theta dtheta J2 psi dpsi L MD mg g  'real'


%% Symbolische Lösung Lagrange-Funktion

% Verallgemeinerte Koordinaten und Ableitungen
q  = [theta,psi];
dq = [dtheta,dpsi];

% kinetische Energien
T_kin = J1*dtheta^2/2 + J2*(dtheta+dpsi)^2/2 - ...
        L*MD*dtheta*(dtheta+dpsi)*cos(theta);

% potentielle Energien
U_pot = +mg*L*g*cos(theta) - MD*g*cos(theta+psi);
U_pot = -mg*L*g*cos(theta) + MD*g*cos(theta+psi);



% Lagrange-Funktion
L = T_kin - U_pot;

%% Ausgabe der Ergebnisse
fprintf("------------------------------------\nLagrange-Funktion L:");
fprintf("\n------------------------------------\n");
fprintf(1, 'T = %s \n', T_kin);
fprintf(1, 'U = %s \n', U_pot);
fprintf(1, 'L = T - U = %s \n', L);
fprintf("\n");

% Externe Kräfte 
Q = [0,0];

% Berechnet Euler-Lagrange-Gleichungen:
EQ = EulerLagrange(q,dq,L,Q,true);



%% Numerische Lösung

%Erwachsener
L       = 4.0;              % Länge Schaukel
hK      = 0.55;             % Länge Oberkörper
hB      = 0.75;             % Länge Beine
mg      = 100;              % Gesamtmasse
mK      = 0.40*mg;          % Masse Oberkörper
mB      = 0.30*mg;          % Masse Beine
%Kind
% L       = 3.0;            % Länge Schaukel
% hK      = 0.30;           % Länge Oberkörper
% hB      = 0.40;           % Länge Beine
% mg      =  50;            % Gesamtmasse
% mK      = 0.30*mg;        % Masse Oberkörper
% mB      = 0.30*mg;        % Masse Beine

mR      = mg-mK-mB;         % Masse Rumpf
J1      = mg*L^2;           % Trägheitsmoment Rumpf
J2      = mB*hB^2+mK*hK^2;  % Trägheitsmoment Kopf + Beine
Jg      = J1+J2;            % Gesamtträgheitsmoment
MD      = mB*hB-mK*hK;      % Drehmoment bzgl. Sitzfläche
g       = 9.81;             % Erdbeschleunigung

psi0  = deg2rad(45);        % Anfangswinkel Körperneigung
K0    = mg*L*g-MD*g*(1-0.25*psi0^2);
J0    = J1+J2-2*L*MD*(1-0.25*psi0^2);
omeg0 = sqrt(K0/J0);        % Eigenfrequenz
% omeg0 = 2;                % Eigenfrequenz
omegA = 1*omeg0;            % Frequenz der Pendellängenänderung

F     = psi0*(omegA^2*J2+MD*(g-omeg0^2*L)*(1-0.125*psi0^2))/J0;
A     = -MD*g*psi0^2/4/J0;
B     = L*MD*omegA*psi0^2/J0;
C     = L*MD*psi0^2/2/J0;
dpsi0   = 0.0;              % Anfangsgeschwindigkeit Körperneigung
phi0    = asin(-dpsi0/psi0/omegA);

TA    = 2*pi/omeg0;         % Schwingungsperiode (in sec)
tmax  = 120;

fprintf('\n ');
fprintf('\n omega_0 = %4.2f', omeg0);
fprintf('\n omega_A = %4.2f', omegA);
fprintf('\n ');
fprintf('\n J1 = %8.3f', J1);
fprintf('\n J2 = %8.3f', J2);
fprintf('\n J0 = %8.3f', J0);
fprintf('\n K0 = %8.3f', K0);
fprintf('\n ');
fprintf('\n F  = %+5.4f', F);
fprintf('\n A  = %+5.4f', A);
fprintf('\n B  = %+5.4f', B);
fprintf('\n C  = %+5.4f', C);



%% Berechnungen

% Anfangswerte
theta0    = 0.0;       % Null
theta0dot = 0.0;
AB=[theta0;theta0dot]; % AB für ode45

% Numerische Lösung Volle Lagrangegleichung
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-7);
% Numerische Lösung LGL
[t1,x1]=ode45(@dgl_schaukel,[0.0,tmax],AB,opt,J1,J2,mg,L,g,MD,omegA,psi0,phi0); 
minx1=max(rad2deg(min(x1(:,1))),-180);
maxx1=min(rad2deg(max(x1(:,1))),180);

% Numerische Lösung linearisierte LGL
[t2,x2]=ode45(@dgl_schaukel1,[0.0,tmax],AB,opt,J1,J2,mg,L,g,MD,omegA,psi0,phi0); 
minx2=max(rad2deg(min(x2(:,1))),-180);
maxx2=min(rad2deg(max(x2(:,1))),180);

% Näherungslösung Reihenentwicklung
[t3,x3]=ode45(@dgl_schaukel2, [0.0,tmax],AB,opt, omeg0, omegA, F, A, B, C); 
minx3=rad2deg(min(x3(:,1)));
maxx3=rad2deg(max(x3(:,1)));

%% Graphische Ausgabe

% Volle Lagrangegleichung
figure();
subplot(3,1,1)
plot(t1,wrapTo180(rad2deg(x1(:,1))),'Color',Colors(2,:), 'LineWidth',2);
hold on
axis([0 tmax minx3 maxx3])
grid on
ylabel('$\theta$ in Grad','interpreter','latex','FontSize',14)
xlabel('$t$ in s','interpreter','latex','FontSize',14)
h=title(' Numerisch Lsg. Lagrange Gl.');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

% Numerische Lösung linearisierte LGL
subplot(3,1,2)
plot(t2,wrapTo180(rad2deg(x2(:,1))),'Color',Colors(3,:), 'LineWidth',2);
hold on
axis([0 tmax minx3 maxx3])
xlabel('$t$ in s','interpreter','latex','FontSize',14)
ylabel('$\theta$ in Grad','interpreter','latex','FontSize',14)
grid on
h=title('Numerisch Lsg. linearisierte LGL');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

% Näherungslösung Reihenentwicklung
subplot(3,1,3)
plot(t3,rad2deg(x3(:,1)),'Color',Colors(4,:), 'LineWidth',2);
hold on
axis([0 tmax minx3 maxx3])
xlabel('$t$ in s','interpreter','latex','FontSize',14)
ylabel('$\theta$ in Grad','interpreter','latex','FontSize',14)
grid on
h=title('Näherung Reihenentwicklung');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

figure()
hp(1)=plot(t3,rad2deg(x3(:,1)),'Color',Colors(4,:), 'LineWidth',2);
hold on
hp(2)=plot(t2,rad2deg(x2(:,1)),'Color',Colors(3,:), 'LineWidth',2);
hp(3)=plot(t1,rad2deg(x1(:,1)),'Color',Colors(2,:), 'LineWidth',2);
if tmax <= 50
    line([0 tmax],[0 maxx3],'Color',Colors(4,:), 'LineWidth',2);
    line([0 tmax],[0 maxx2],'Color',Colors(3,:), 'LineWidth',2);
    line([0 tmax],[0 maxx1],'Color',Colors(2,:), 'LineWidth',2);
end
xlabel('$t$ in s','interpreter','latex','FontSize',14)
ylabel('$\theta$ in Grad','interpreter','latex','FontSize',14)
grid on
hlgd=legend(hp,'Reihenentwicklung','Lsg. linear. LGL','Lsg. Lagrange GL',...
          'location','best');
set(hlgd,'FontSize',14,'FontWeight','normal'); 
legend box off
h=title('Vergleich der Lösungen');
set(gca,'FontSize',16);
set(h,'FontSize',14,'FontWeight','normal'); 
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktioen

% Volle Lagrangegleichung
function dY = dgl_schaukel(t, Y, J1,J2,mg,L,g,MD,omegA,psi0,phi0)
    % Y(1)-Winkel theta
    % Y(2)-Winkelgeschwindigkeit thetadot
    psi   =  psi0*cos(omegA*t+phi0);
    dpsi  = -psi0*omegA*sin(omegA*t+phi0);
    ddpsi = -psi*omegA^2;
    dY    = [Y(2);...
        ((- mg*L*g  + MD*g*cos(psi))*sin(Y(1))...
         - J2*ddpsi + MD*L*cos(psi)*ddpsi ...
         - MD*L*sin(psi)*dpsi^2 + MD*g*sin(psi)*cos(Y(1)) ...
         - 2*L*MD*sin(psi)*dpsi*Y(2)) ...
         /(J1+J2-2*MD*L*cos(psi))]; 
end

% Linearisierte Bewegungsgleichung
function dY = dgl_schaukel1(t, Y, J1,J2,mg,L,g,MD,omegA,psi0,phi0)
    % Y(1)-Winkel theta
    %Y(2)-Winkelgeschwindigkeit thetadot
    psi   =  psi0*cos(omegA*t+phi0);
    dpsi  = -psi0*omegA*sin(omegA*t+phi0);
    ddpsi = -psi*omegA^2;
    dY    = [Y(2);...
        ((- mg*L*g  + MD*g*cos(psi))*Y(1)...
         - J2*ddpsi + MD*L*cos(psi)*ddpsi ...
         - MD*L*sin(psi)*dpsi^2 + MD*g*sin(psi)*1 ...
         - 2*L*MD*sin(psi)*dpsi*Y(2)) ...
         /(J1+J2-2*MD*L*cos(psi))]; 
end


% Reihenentwicklung
function dY = dgl_schaukel2(t, Y, omeg0, omegA, F, A, B, C)
    % Y(1)-Winkel theta
    % Y(2)-Winkelgeschwindigkeit thetadot
    arg   =  omegA*t;
    dY    = [Y(2);...
        (- omeg0^2*Y(1)     + F*cos(arg) ...
         + A*cos(2*arg)*Y(1) + B*sin(2*arg) *Y(2))...
         /(1-C*cos(2*arg))]; 
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


