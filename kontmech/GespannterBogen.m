% -------------------------------------------------------------------------
% GespannterBogen.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Schwingungsfrequenz eines gespannten Bogens
% in Abhängigkeit von der anfänglichen Biegung.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Materialparameter
M = 0.118;         % Masse des Bogens in kg
EI = 3.09;         % E*I in kg*m^2
L = 0.52;          % Länge des Bogens in m
k = 2.1e4;         % Federkonstante der Saite in N/m
T1 = 100;          % Kraft der Saite, Ruhelage des gesp. Bogens in N

% Diskretisierung des Bogens
Ndim = 4;              % Dimension des DGL-Systems
Smin = 0.;             % Start der Bogenlänge
Smax = L;              % Bogenlänge in m
Nint = 500;            % Anzahl der Integrationsschritte

Ndim = 4;              % Dimension des DGL-Systems
Smin = 0.;             % Start der Bogenlänge
Smax = L;              % Bogenlänge in m
Nint = 2000;           % Anzahl der Integrationsschritte
ds   = Smax/(Nint-1);  % Schrittweite der Integration

%% Nullstellensuche: Finde die passenden Frequenzen

% Anfangswerte fuer die Suche ohne Saite
omega = 2*pi*35.;
dnds = -6;
n0 = 1.;

% Angriffspunkte der Saite (Abstand vom Rand) auf nicht erreichbare
% Werte gesetzt -> Es gibt keine Saite
N1 = -5; N2 = -5;

% Suche die Frequenz ohne Saite
P0(1) = omega;
P0(2) = dnds;
theta(1) = 0.;
SFUNC = @(P)FUNC(P,M,EI,L,k,theta(1),T1,ds,n0,Nint,Ndim,N1,N2);
P = fsolve(SFUNC,P0);
omega = P(1);
dnds = P(2);

% Speichere den Wert
omegaO = omega;

%% Nun mit Saite und Variation des Winkes theta

% Angriffspunkte der Saite
Ndist = 10;
N1 = 1+Ndist;
N2 = Nint-Ndist;

theta1 = 0.;           % Beginn des Intervalls der anfänglichen Biegung
theta2 = 25/180.*pi;   % Ende des Intervalls 
Ntheta = 100;          % Anzahl der betrachteten Biegungen

% T1 = 20
T1 = 20;               % Kraft der Saite, Ruhelage des gesp. Bogens in N
for j=1:Ntheta+1
    theta(j) = theta1+(theta2-theta1)/Ntheta*(j-1);

    P0(1) = omega;
    P0(2) = dnds;
    SFUNC = @(P)FUNC(P,M,EI,L,k,theta(j),T1,ds,n0,Nint,Ndim,N1,N2);
    P = fsolve(SFUNC,P0);
    omega = P(1);
    dnds = P(2);
    
    % Speichere den Wert
    omega20(j) = omega;
end

% Wiederholung für T1 = 50
omega = 48*2*pi;
dnds = -0.6;
T1 = 50;         
for j=1:Ntheta+1
    theta(j) = theta1+(theta2-theta1)/Ntheta*(j-1);

    P0(1) = omega;
    P0(2) = dnds;
    SFUNC = @(P)FUNC(P,M,EI,L,k,theta(j),T1,ds,n0,Nint,Ndim,N1,N2);
    P = fsolve(SFUNC,P0);
    omega = P(1);
    dnds = P(2);
    
    % Speichere den Wert
    omega50(j) = omega;
end

% Wiederholung für T1 = 100
omega = 48*2*pi;
dnds = -0.6;
T1 = 100;         
for j=1:Ntheta+1
    theta(j) = theta1+(theta2-theta1)/Ntheta*(j-1);

    P0(1) = omega;
    P0(2) = dnds;
    SFUNC = @(P)FUNC(P,M,EI,L,k,theta(j),T1,ds,n0,Nint,Ndim,N1,N2);
    P = fsolve(SFUNC,P0);
    omega = P(1);
    dnds = P(2);
    
    % Speichere den Wert
    omega100(j) = omega;
end

%% Frequenzverlauf in einer Abbildung
figure
plot(theta(:)/pi*180,omega20(:)/2./pi, 'color', Colors(2,:), ...
    'LineWidth',2)
hold on
plot(theta(:)/pi*180,omega50(:)/2./pi, 'color', Colors(3,:), ...
    'LineWidth',2)
hold on
plot(theta(:)/pi*180,omega100(:)/2./pi, 'color', Colors(4,:), ...
    'LineWidth',2)
hold on
plot(theta(:)/pi*180,ones(size(theta))*omegaO/2./pi, ...
    'color', Colors(1,:),'LineWidth',2)
ylabel('\itf')
xlabel('\theta')
grid on
legend('{\itf}_{20}', '{\itf}_{50}', '{\itf}_{100}', '{\itf}_{0}', ...
    'location','northwest','numcolumns',1);
legend box off
set(gca,'FontSize',16,'FontName','Times');


%% Berechnung der Auslenkung für den letzten Parameterwert
Afrei = M/Nint*omega^2/(EI*L/Nint);
Bfrei = 0;
ASaite = (M/Nint*omega^2-2*k*sin(theta(Ntheta+1))^2)/(EI*L/Nint);
BSaite = T1/(EI*L/Nint)*cos(theta(Ntheta+1));

Y0(1) = n0;
Y0(2) = dnds;
Y0(3) = 0.;
Y0(4) = 0.;

% Parameter setzen
h = ds; 
X = (0:Nint-1)*h;          % Intervall in der unabhängigen Variable
Y = zeros(Ndim,length(X)); % abhängige Variable(n)
Y(:,1) = Y0(:);            % Setzen der Anfangswerte

%Integration
for i = 1:Nint-1
  if (i==N1 | i==N2)
      A = ASaite;
      B = BSaite;
  else
      A = Afrei;
      B = Bfrei;
  end
  dY =@(X,Y)DGL(X,Y,A,B,Ndim);
  k1 = dY(X(i),Y(:,i));
  k2 = dY(X(i)+.5*h,Y(:,i)+.5*k1(:)*h);
  k3 = dY(X(i)+.5*h,Y(:,i)+.5*k2(:)*h);
  k4 = dY(X(i)+h,Y(:,i)+k3(:)*h);
  Y(:,i+1) = Y(:,i) + (k1(:)+2*k2(:)+2*k3(:)+k4(:))*h/6.;
end

%% Abbildungen der Auslenkung für den letzten Parameterwert

% n, dn/ds
figure
plot(X*100,Y(1,:), 'color', Colors(1,:),'LineWidth',2)
hold on;
plot(X*100,Y(2,:), 'color', Colors(2,:),'LineWidth',2)
hold off;
ylabel('{\itn}, d{\itn}/d{\its}')
xlabel('{\its} in cm','FontSize',14)
grid on
legend('\itn', 'd{\itn}/d{\its}','location','southeast','numcolumns',1);
legend box off
set(gca,'FontSize',16,'FontName','Times');

% d^2n/ds^2, d^3n/ds^3
figure
plot(X,Y(3,:), 'color', Colors(3,:),'LineWidth',2)
hold on;
plot(X,Y(4,:), 'color', Colors(4,:),'LineWidth',2)
hold off;
ylabel('d^2{\itn}/d{\its}^2, d^3{\itn}/d{\its}^3')
xlabel('{\it s} in cm','FontSize',14)
grid on
legend('d^2{\itn}/d{\its}^2', 'd^3{\itn}/d{\its}^3', ...
    'location','southeast','numcolumns',1);
legend box off
set(gca,'FontSize',16,'FontName','Times');

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%%
% Funktionen

% Integration für einen Satz an Anfangswerten
function F = FUNC(XF,M,EI,L,k,theta,T1,ds,n0,Nint,Ndim,N1,N2)

  % Extrahieren der Eingabewerte
  omega = XF(1);
  Y0(2) = XF(2);

  % Setzen der restlichen Angangsbedingungen
  Y0(1) = n0;
  Y0(3) = 0.;
  Y0(4) = 0.;

  % Konstanten für die Differentialgleichung
  Afrei = M/Nint*omega^2/(EI*L/Nint);
  Bfrei = 0;
  ASaite = (M/Nint*omega^2-2*k*sin(theta)^2)/(EI*L/Nint);
  BSaite = T1/(EI*L/Nint)*cos(theta);

  % Integration
  h = ds; 
  X = (0:Nint-1)*h;          % Intervall in der unabhängigen Variable
  Y = zeros(Ndim,length(X)); % abhängige Variable(n)
  Y(:,1) = Y0(:);            % Setzen der Anfangswerte
  for i = 1:Nint-1
    if (i==N1 | i==N2)
        A = ASaite;
        B = BSaite;
    else
        A = Afrei;
        B = Bfrei;
    end
    dY =@(X,Y)DGL(X,Y,A,B,Ndim);
    k1 = dY(X(i),Y(:,i));
    k2 = dY(X(i)+.5*h,Y(:,i)+.5*k1(:)*h);
    k3 = dY(X(i)+.5*h,Y(:,i)+.5*k2(:)*h);
    k4 = dY(X(i)+h,Y(:,i)+k3(:)*h);
    Y(:,i+1) = Y(:,i) + (k1(:)+2*k2(:)+2*k3(:)+k4(:))*h/6.;
  end

  F(1) = Y(3,Nint);
  F(2) = Y(4,Nint);
end


% Differentialgleichungssystem
function dY = DGL(X,Y,A,B,N)
dY  =  zeros(N,1); % richtige Dimension von dY

dY(1)  = Y(2);
dY(2)  = Y(3); 
dY(3)  = Y(4);
dY(4)  = A*Y(1)+B*abs(Y(2)); 

end
% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
