% -------------------------------------------------------------------------
% FallenderStab02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Dynamik des Fallender Stab auf Basis der
% Lagrange-Gleichungen  
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

L       =  0.254;           % Länge Stab in  m
H       =  L/2;             % Schwerpunktlage bei homogenen Stab  m
                            % kann auch verändert werden für inhomogenen S.
m       =  0.0554;          % Masse des Stabs in kg
mu1     =  0.10;            % Reibungskoeffizient 1 
mu2     =  0.30;            % Reibungskoeffizient 2 
mu3     =  0.60;            % Reibungskoeffizient 3 
JS      =  m*L^2/12;        % Trägheitsmoment Schwerpunkt homogener Stab
JP      =  m*L^2/3;         % Trägheitsmoment Fusspunkt homogener Stab
                            % kann auch verändert werden für inhomogenen S.
g       =  9.81;            % g in m/s^2 
theta0  =  deg2rad(1);      % Anfangswinkel 
dtheta0 =  0;               % Anfangsgeschwindigkeit rad/s
tmax    =  2.0;             % Fallzeit in s
tspan   =  linspace(0.0,tmax,200);

fprintf('\n ');
fprintf('\n L   = %8.2f m', L);
fprintf('\n m   = %8.2f kg', m);
fprintf('\n JS  = %8.2f g cm^2', JS*10^7);
fprintf('\n JP  = %8.2f g cm^2', JP*10^7);
fprintf('\n mu1 = %8.2f ', mu1);
fprintf('\n mu2 = %8.2f ', mu2);
fprintf('\n mu3 = %8.2f ', mu3);
fprintf('\n g   = %8.2f m/s^2', g);
fprintf('\n tmax= %8.2f m/s', tmax);
fprintf('\n ');


%% Berechnungen

% Anfangswerte
AB=[theta0;dtheta0]; % AB für ode45
mu  = [mu1, mu2, mu3];
omega02=m*g*H/JP;

% Berechnung der Slip-Winkel

theta = linspace(deg2rad(1),deg2rad(90),90);
FH   = m*H*omega02*sin(theta).*(3*cos(theta)-2*cos(theta0));
FN   = m*g-m*H*omega02*(1-3*cos(theta).*cos(theta)+2*cos(theta).*cos(theta0));
% Bestimmung Slip-Winkel thetaS
myfun = @(x,c1,c2) (sin(x).*(9*cos(x)-6*cos(c2))./...
                   (1+9*cos(x).*cos(x)-6*cos(x).*cos(c2))) - c1;  
                 % parameterized function               
thetaC = acos(2*cos(theta0)/3);                
for k=1:3
    fun = @(x) myfun(x,mu(k),theta0);    % function of x alone
    thetaS(k) = fzero(fun,theta0);
    if (thetaS(k) < 0) || (thetaS(k) > thetaC)
        thetaS(k) = thetaC;
        P2.forward(k) = true;
    else 
        P2.forward(k) = false;
    end        
end


%%
% Dynamikberechnung für homogenen Stab

% Parameter für DGL Phase1
P1.om     = 3*g/(2*L);
P1.thetaS = thetaS; 

% Parameter für DGL Phase2
P2.mu  = mu1;
P2.index  = 1;
P2.g  = g;
P2.H  = H;


t1 = NaN(3,100);
z1 = NaN(3,100);
for k=1:3
    P2.mu = mu(k);
    P2.index = k;
    P1.index = k;
    % Numerische Lösung LGL
    mustr(k,:) = strcat('\mu = ',string(num2str(P2.mu,'% 5.2f')));
    options = odeset('AbsTol',1.e-7,'RelTol',1.e-5,'events',@MyEvent1);
    [t, Y, TE, YE, IE] = ode45(@dgl_FallingRod1,tspan, AB,options,P1); 
    switch k
    case 1
        t1   = t;
        TE1  = TE;
        z1   = rad2deg(Y(:,1));
        xF1  = zeros(length(t1),1);
    case 2
        t2   = t;
        TE2  = TE;
        z2   = rad2deg(Y(:,1));
        xF2  = zeros(length(t2),1);
    case 3
        t3   = t;
        TE3  = TE;
        z3   = rad2deg(Y(:,1));
        xF3  = zeros(length(t3),1);
    end
    AB2 =[YE(1); YE(2);0;0]; % AB für 2. Phase ode45
    options = odeset('AbsTol',1.e-7,'RelTol',1.e-5,'events',@MyEvent2);
    [t, Y] = ode45(@dgl_FallingRod2,tspan, AB2,options,P2); 
    switch k
        case 1
            t12   = TE1+t;
            z12   = rad2deg(Y(:,1));
            xP1 = Y(:,3);
        case 2
            t22   = TE2+t;
            z22   = rad2deg(Y(:,1));
            xP2 = Y(:,3); 
        case 3
            t32   = TE3+ t;
            z32   = rad2deg(Y(:,1));
            xP3 = Y(:,3); 
     end
end

%% 
% Graphische Ausgabe

tplot = 1.2*max([t1;t2;t3]);
fig=figure();
% set(fig,'defaultAxesColorOrder',[Colors(3,:); Colors(2,:)]);

% Fallwinkel
subplot(1,2,1)
h=title('Numerische Lösung Lagrange-Gl.');
hold on
lp(1) = plot(t1,abs(z1));
lp(2) = plot(t2,abs(z2));
lp(3) = plot(t3,abs(z3));
for k=1:3 
    set(lp(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(k));
end

temp(1) = plot(t12,abs(z12));
temp(2) = plot(t22,abs(z22));
temp(3) = plot(t32,abs(z32));
for k=1:3 
    set(temp(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(k));
end
axis([0 tplot 0 90]);
grid on
ylabel('Winkel \theta  in °','FontSize',14)
xlabel('t in s','FontSize',14)
legend(lp, mustr, 'location','northwest');
legend box off
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

% Verschiebung Fußpunkt
subplot(1,2,2)
hold on;
rp(1) = plot(t1, xF1*1000, 'Color',Colors(2,:), 'LineWidth',2,'LineStyle',Style(1));
rp(2) = plot(t2, xF2*1000, 'Color',Colors(3,:), 'LineWidth',2,'LineStyle',Style(2));
rp(3) = plot(t3, xF3*1000, 'Color',Colors(4,:), 'LineWidth',2,'LineStyle',Style(3));
rp(1) = plot(t12, xP1*1000, 'Color',Colors(2,:), 'LineWidth',2,'LineStyle',Style(1));
rp(2) = plot(t22, xP2*1000, 'Color',Colors(3,:), 'LineWidth',2,'LineStyle',Style(2));
rp(3) = plot(t32, xP3*1000, 'Color',Colors(4,:), 'LineWidth',2,'LineStyle',Style(3));
grid on
axis([0 tplot 1.2*min(xP1*1000) 1.2*max(xP3*1000)]);
xlabel('t in s','FontSize',14)
ylabel('Fußpunktverschiebung x_P in mm','FontSize',14)
legend(rp, mustr, 'location','northwest');
legend box off
grid on
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

function dY = dgl_FallingRod1(t,Y,P1)
    % q- Winkel theta(t), 
    % dq- Winkelgeschwindigkeit dtheta(t)
    dY    = [Y(2);...    
             P1.om*sin(Y(1))];
end

function dY = dgl_FallingRod2(t,Y,P2)  % Rechnung für homogenen Stab
    % q1 - Winkel theta(t), 
    % dq1- Winkelgeschwindigkeit dtheta(t)
    % q2 - Fusspunkt xp(t), 
    % dq2- Geschwindigkeit dxp(t)
    q1  = Y(1);
    dq1 = Y(2);
    q2  = Y(3);
    dq2 = Y(4);
    if P2.forward
        muR = -P2.mu;
    else
        muR = P2.mu;
    end
    b = 1./(3*sin(q1)-3*muR*cos(q1));
    alpha = (P2.g/P2.H - dq1.*dq1.*cos(q1))./(sin(q1)+b);
    dY    = [dq1;...    
            alpha;...
            dq2;...
            dq1.*dq1*P2.H*sin(q1)-alpha*P2.H.*(cos(q1)-b*muR)];
end

 
% Ereignisfunktion

function [value,isterminal,direction] = MyEvent1(t,Y,P1)
    % Ereignisfunktion bis Eintreten des Gleitens bei Slip Winkel thetaS bzw.
    % thetaC
    %     FH   = (3/4)*P1.m*P1.g*sin(Y(1)).*(3*cos(Y(1))-2*cos(P1.wi));
    %     FN   = (1/4)*P1.m*P1.g*(1+9*cos(Y(1)).*cos(Y(1))-6*cos(Y(1))*cos(P1.wi));
    value = (P1.thetaS(P1.index)-Y(1));
    % erkenne Slip Winkel
    isterminal = 1;   % stop Integration
    direction = 0;    % negative Richtung
end

function [value,isterminal,direction] = MyEvent2(t,Y,P1)
    value = pi/2-Y(1); % erkenne Winkel 90°
    isterminal = 1;    % stop Integration
    direction = 0;     % beliebige Richtung
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------



