% -------------------------------------------------------------------------
% Wilberforce01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Wilberforce Pendel
% Berechnung der Dynamik des Wilberforce-Pendels auf 
% Basis von Lagrange-Gleichungen (Resonanzfall)
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Numerische Lösung
% Parameter
mF = 0.5164;            % Masse Federpendel, kg
kF = 2.8;               % Federkonstante, N/m
LF = 1.765;             % Länge Feder, m
JT = 1.45e-4;           % Trägheitsmoment der Torsionsstange, kg m2
DT = 7.9e-4;            % Torsionsmomemt, Nm/rad
kappa  = 9.27e-3;       % Kopplungskonstante 1, N
omegaF = sqrt(kF/mF)    % Eigenfrequenz Feder
omegaT = sqrt(DT/JT)    % Eigenfrequenz Torsion
omegaB = kappa/2/omegaT/sqrt(mF*JT)    % Schwebungsfrequenz 
TF = 2*pi/omegaF
TB = 2*pi/omegaB
etaF = 0;               % Reibungsterm Feder
etaT = 0;               % Reibungsterm Torsion

A = 0.03;               % Amplitude Antrieb Feder
B = 0;                  % Amplitude Antrieb Torsion
omega = 1.0*omegaF;     % Antriebsfrequenz

epsilonF=0.00;          % nichtlinearer Anteil Feder
epsilonT=0.0;           % nichtlinearer Anteil Torsion


%% 
% Anfangswerte
Y10 = 0.3;              % z     Federposition
Y20 = 0.0;              % dz    Federgeschwindigkeit
Y30 = 0.0;              % phi   Torsionswinkel
Y40 = 0;                % dphi  Torsionswinkelgeschwindigkeit
Y0 = [Y10 Y20 Y30 Y40]; % Anfangsvektor


% lineare DGL dz = AA*z, exakte lineare Lösung
% AA = [0             1         0           0;
%      -omegaF^2     -etaT/mF  -kappa/2     0;
%       0             0         0           1;
%      -kappa/2       0        -omegaT^2   -etaF/JT]
% eig_AA = eig(AA)
% [Y,X] = initial(AA,zeros(4,1),eye(4),zeros(4,1),Y0,t);
% figure('Exakte Lösung
% subplot(211),plotyy(t,Y(:,1),t,Y(:,3)), ylabel('z, phi'), 
% title('Lineare LÃ¶sung, exakt')
% subplot(212),plotyy(t,Y(:,2),t,Y(:,4)), ylabel('dz, dphi')

%% 
% 
%Numerische Lösung volle Lagrangegleichung Linear
%
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-8);

tmax = 60;             % Simulationszeit
tspan = linspace(0,tmax,600);

% Numerische Lösung LGL
[t,Y]=ode45(@dgl_Wilberforce,tspan,Y0,opt,mF,omegaF,etaF,kappa,etaT,...
                             omegaT,JT,0,0,omega,0,0); 

z = Y(:,1);
dz= Y(:,2);
phi  = Y(:,3);
dphi = Y(:,4);

% Graphische Ausgabe
figure(1)
subplot(211)
yyaxis left
plot(tspan,z)
ylabel ('\it z \rm in m')
yyaxis right
plot(tspan,phi/pi)
ylabel ('\phi in \pi')
h=title('Lineare Lösung mit ODE45');
% subplot(212),plotyy(t,dz,t,dphi), ylabel('dz, dphi')
grid on
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

ET= 0.5*JT*dphi.*dphi+0.5*DT*phi.*phi;
EF= 0.5*mF*dz.*dz+0.5*kF*z.*z;
subplot(212),plot(t,EF,t,ET,t,ET+EF+0.5*kappa*z.*phi)
xlabel('\it t \rm in s'), ylabel('\it E \rm in Nm'),
h=title('Lineare Lösung, Energien');  
legend('E_F','E_T','E_{ges}','location', 'west')
legend box off
grid on
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

%%
% Analytische Lösung
%
tmax = 60;                  % Simulationszeit
ta = linspace(0,tmax,600);
z0   = 0.3;                 % z     Federposition
phi0 = 0.0;                 % phi   Torsionswinkel

omega1 = 0.5*(omegaF^2+omegaT^2+sqrt((omegaF^2-omegaT^2)^2+kappa^2/mF/JT));
omega1 = sqrt(omega1)

omega2 = 0.5*(omegaF^2+omegaT^2-sqrt((omegaF^2-omegaT^2)^2+kappa^2/mF/JT));
omega2 = sqrt(omega2)

diff12  = (omega1^2-omega2^2);
diff1T  = (omega1^2-omegaT^2);
diff2T  = (omega2^2-omegaT^2);

if kappa >0
    Z   = (z0/diff12)*(diff1T*cos(omega1*ta)-diff2T*cos(omega2*ta))...
		- (2*JT*phi0/kappa/diff12)*diff1T*diff2T...
        *(cos(omega1*ta)-cos(omega2*ta));
else
    Z   = (z0/diff12)*(diff1T*cos(omega1*ta)-diff2T*cos(omega2*ta));
end
Phi = (kappa*z0/2/JT/diff12)*(cos(omega1*ta)-cos(omega2*ta));

figure(2)
subplot(211)
yyaxis left
plot(t,Z(:)-z(:))
ylabel ('\it Differenz z \rm in m')
yyaxis right
plot(ta,(Phi(:)-phi(:))/pi)
ylabel ('Diffrenz \phi in \pi')
h=title('Abweichung Numerischer von analytischer Lösung');
grid on
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

% Anfangsbedingung so, dass nur eine Normalschwingung angeregt wird
z0     = 0.005;                     % z     Federposition
phi0   = sqrt(mF/JT)*z0 ;          % phi   Torsionswinkel

if kappa >0
    Z   = (z0/diff12)*(diff1T*cos(omega1*ta)-diff2T*cos(omega2*ta))...
		- (2*JT*phi0/kappa/diff12)*diff1T*diff2T...
        *(cos(omega1*ta)-cos(omega2*ta));
else
    Z   = (z0/diff12)*(diff1T*cos(omega1*ta)-diff2T*cos(omega2*ta));
end
Phi = (kappa*z0/2/JT/diff12)*(cos(omega1*ta)-cos(omega2*ta));
vecn1 = [1 2*JT*diff1T/kappa];
vecn2 = [1 2*JT*diff2T/kappa];
N1  = (1/diff12)*(kappa*z0/2/JT - diff2T*phi0);
N2 = -(1/diff12)*(kappa*z0/2/JT - diff1T*phi0);
XPhi    = N1*cos(omega1*ta) + N2*cos(omega2*ta);
XZ      = (N1*cos(omega1*ta)*diff1T + N2*cos(omega2*ta)*diff2T)*2*JT/kappa;
X(1,:)   = N1*cos(omega1*ta).*vecn1(1) + N1*cos(omega1*ta).*vecn2(1);
X(2,:)   = N1*cos(omega1*ta).*vecn1(2) + N1*cos(omega1*ta).*vecn2(2);

subplot(212)
yyaxis left
plot(ta,XZ)
ylabel ('\it z \rm in m')
yyaxis right
plot(ta,XPhi)
xlabel('\it t \rm in s')
ylabel ('\phi in \pi')
h=title('Normalschwingung');
grid on
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);





%% Numerische Lösung volle LGL mit Nichtlinearität und Anregung
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-6);

Y10 = 0.3;              % z     Federposition
Y20 = 0.0;              % dz    Federgeschwindigkeit
Y30 = 0.0;              % phi   Torsionswinkel
Y40 = 0;                % dphi  Torsionswinkelgeschwindigkeit
Y0 = [Y10 Y20 Y30 Y40]; % Anfangsvektor

etaF = 1e-04;           % Reibungsterm Feder
etaT = 1e-06;           % Reibungsterm Torsion

A = 0.1;                % Amplitude Antrieb Feder
B = 0.0;                % Amplitude Antrieb Torsion
omega = 1.0*omegaF;     % Antriebsfrequenz

epsilonF=0.03;          % nichtlinearer Anteil Feder
epsilonT=0.00;          % nichtlinearer Anteil Torsion

tmax = 60;             % Simulationszeit
t = linspace(0,tmax,600);

% Numerische Lösung LGL
[t,Y]=ode45(@dgl_Wilberforce,t,Y0,opt,mF,omegaF,etaF,kappa,etaT, ...
                             omegaT,JT,A,B,omega,epsilonF,epsilonT); 

% Graphische Ausgabe

z    = Y(:,1);
dz   = Y(:,2);
phi  = Y(:,3);
dphi = Y(:,4);

figure(3)
subplot(211)
yyaxis left
plot(t,z)
ylabel ('\it z \rm in m')
yyaxis right
plot(t,phi/pi)
ylabel ('\phi in \pi')
xlabel('\it t \rm in s')
h=title('Nichtlineare Lösung mit ODE45');
% subplot(212),plotyy(t,dz,t,dphi), ylabel('dz, dphi')
grid on
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);


ET= 0.5*JT*dphi.*dphi+0.5*DT*phi.*phi;
EF= 0.5*mF*dz.*dz+0.5*kF*z.*z;
subplot(212),plot(t,EF,t,ET,t,ET+EF+0.5*kappa*z.*phi)
xlabel('\it t \rm in s'), ylabel('\it E \rm in Nm'),
h=title('Nichtlineare Lösung, Energien');  
legend('E_F','E_T','E_{ges}','location', 'northwest')
legend box off
grid on
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);



% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%% Funktionen
% Lagrange-Gleichung
function dY = dgl_Wilberforce(t,Y,mF,omegaF,etaF,kappa,etaT,omegaT,JT,...
                              A,B,omega,epsilonF,epailonT)
dY = [   Y(2)
        -omegaF^2*(Y(1)+epsilonF*Y(1)^3)-1/2*kappa/mF*Y(3)-etaF/mF*Y(2)... 
              + A*sin(omega*t)
         Y(4)
        -omegaT^2*(Y(3)+epailonT*Y(3)^3)-1/2*kappa/JT*Y(1)-etaT/JT*Y(4)...
              + B*sin(omega*t)
     ];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
