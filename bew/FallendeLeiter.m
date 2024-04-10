% -------------------------------------------------------------------------
% FallendeLeiter.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Dynamik der fallenden Leiter 
% auf Basis Lagrange-Gleichungen mit und ohne Reibung  
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

L       = 5;               % Länge Leiter m
m       = 10;              % Masse in kg
g       = 9.81;            % g in m/s
mu      = 0.3;             % Reibungskoeffizient
tmax    = 5;
% Parameter in LGL
% Anfangswerte so, dass Leiter durch kleine Störungen ins Rutschen kommt
phi0 = deg2rad(87.0);
dphi0= -0.00005;


fprintf('\n ');
fprintf('\n ');
fprintf('\n L    = %8.2f m', L);
fprintf('\n m    = %8.2f kg', m);
fprintf('\n g    = %8.2f m/s', g);
fprintf('\n mu   = %8.2f ', mu);
fprintf('\n phi0 = %8.2f °', rad2deg(phi0));
fprintf('\n dphi0= %8.2f 1/s', dphi0);
fprintf('\n ');



%% Berechnungen ODE45

% Anfangswerte
AB=[phi0;dphi0]; % AB für ode45

P1.m     = m;
P1.g     = g;
P1.L     = L;
P1.mu    = mu;


% Numerische Lösung Lagrangegleichung
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-7,'events',@MyEvent);
% Numerische Lösung LGL ohne Reibung
[t, Y, TE, YE, IE] = ode45(@dgl_Ladder01,[0.0,tmax],AB,opt,P1); 
phi  = Y(:,1);
dphi = Y(:,2); 
kend = length(phi);
alpha0 = dphi0*sqrt(L/g);
FBz   = m*g-m*L*0.5.*(1.5*g.*cos(phi).*cos(phi)/L+sin(phi).*dphi.*dphi);
FWx   = m*L*0.5*cos(phi).*(1.5*g.*sin(phi)/L - dphi.*dphi);
% FBz1  = m*g*0.25*(cos(phi).*cos(phi)+2*sin(phi).*(5*sin(phi) -3*sin(phi0)-alpha0^2));
% FWx1  = m*g*0.25*cos(phi).*(9*sin(phi) - 6*sin(phi0) - 2*alpha0^2);

% Numerische Lösung LGL mit Reibung

if phi0 > atan((1-mu^2)/2/mu)  % Leiter Stabil, kein Rutschen
    fprintf('\n');
    fprintf('\n Reibung mit mu = %s zu groß.',num2str(mu,3));
    fprintf('\n Leiter wird bei den Anfangsbedingungen nicht wegrutschen.');
    fprintf('\n');
    FrictionCase = 0;
else  % Leiter kommt ins Rutschen
    fprintf('\n');
    fprintf('\n Reibung mit mu = %s klein genug.',num2str(mu,3));
    fprintf('\n Leiter wird bei den Anfangsbedingungen wegrutschen.');
    fprintf('\n');
    FrictionCase = 1;
    [t1, Y1, TE, YE, IE] = ode45(@dgl_Ladder02,[0.0,tmax],AB,opt,P1); 
    phi1  = Y1(:,1);
    dphi1 = Y1(:,2); 
    kend1 = length(phi1);   
    FBz1  = (m*g/2*(5-3.*cos(2*phi1)-3*mu.*sign(dphi1).*sin(2*phi1)) - ...
                   m*L.*(2*sin(phi1)-mu.*sign(dphi1).*cos(phi1)).*dphi1.^2)*...
                   (1/2/(2-mu^2));
               
    FWx1  = (m*g/2*(mu*sign(dphi1).*(1-3.*cos(2*dphi1))+3.*sin(2*phi1))-...
             m*L.*(2*cos(phi1)+mu*sign(dphi1).*sin(phi1)).*dphi1.^2)*...
              (1/2/(2-mu^2));
end


%% 
% Graphische Ausgabe

% Zeitentwicklung 
figure();
hold on
p(1) = plot(t, rad2deg(phi),'Color',Colors(3,:), 'LineWidth',2);
if FrictionCase > 0 
    p(2) = plot(t1, rad2deg(phi1),'Color',Colors(4,:), 'LineWidth',2,...
    'LineWidth',2,'LineStyle',Style(2));
    axis([0 1.2*max(t1) 0 1.2*max(rad2deg(phi))]);

else
    p(2) = line([0 t(kend)],[rad2deg(phi0),rad2deg(phi0)],'Color',...
        Colors(4,:),'LineWidth',2,'LineStyle',Style(2));   
    axis([0 1.2*max(t) 0 1.2*max(rad2deg(phi))]);
end
grid on
ylabel('\phi °','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('\phi (t) ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p,'ohne Reibung','mit Reibung', 'location','northeast',...
       'NumColumns',2);
legend box off
set(gca,'FontSize',16);

[maxFWx,Imax] = max(FWx);
Izero = find(abs(FWx/m/g) < 0.01);
if Izero > 0
    phizero = rad2deg(phi(Izero(1)));
    phimax = rad2deg(phi(Imax));
else
    phizero = rad2deg(phi(1));
    phimax = rad2deg(phi(Imax));
    Izero  = 0;
end

if FrictionCase > 0
    [maxFWx1,Imax1] = max(FWx1);
    Izero1 = find(abs(FWx1/m/g) < 0.02);
    if Izero1 > 0
        phizero1 = rad2deg(phi1(Izero1(1)));
        phimax1  = rad2deg(phi1(Imax1));
    else
        phizero1 = rad2deg(phi1(1));
        phimax1  = rad2deg(phi1(Imax1));
        Izero1  = 0;
    end
else
    Izero1  = 0;
    phizero1 = rad2deg(phi0);
    phimax1  = rad2deg(phi0);

 end
    

% Kräfte 
figure()
p1(1) = plot(rad2deg(phi),FBz/m/g,'Color',Colors(2,:), 'LineWidth',2);
hold on
p1(2) = plot(rad2deg(phi),FWx/m/g,'Color',Colors(3,:), 'LineWidth',2);
if Izero >0 
    p1(3) = line([phizero,phizero],[-2 0],'Color',Colors(3,:),...
    'LineWidth',1,'LineStyle',Style(1));
    p1(4) = line([phimax,phimax],[-2 maxFWx/m/g],'Color',Colors(5,:),...
        'LineWidth',1,'LineStyle',Style(1));
else
    p1(3) = line([0,0],[0 0],'Color',Colors(3,:),...
    'LineWidth',1,'LineStyle',Style(1));   
    p1(4) = line([phimax,phimax],[-2 maxFWx/m/g],'Color',Colors(5,:),...
        'LineWidth',1,'LineStyle',Style(1));
end

if FrictionCase >0
    p1(5) = plot(rad2deg(phi1), FBz1/m/g,'Color',Colors(2,:),...
        'LineWidth',2,'LineStyle',Style(2));   
    p1(6) = plot(rad2deg(phi1), FWx1/m/g,'Color',Colors(3,:),...
        'LineWidth',2,'LineStyle',Style(2));   
else
    p1(5) = line([rad2deg(phi0),rad2deg(phi0)],[-2 2],'Color',Colors(5,:),...
        'LineWidth',1,'LineStyle',Style(2));
end
axis([0 min(90,rad2deg(phi0)+10) -2 2]);
set (gca,'Xdir','reverse')
if Izero1*FrictionCase >0
    p1(7) = line([phizero1,phizero1],[-2 0],'Color',Colors(3,:),...
    'LineWidth',1,'LineStyle',Style(2));
    p1(8) = line([phimax1,phimax1],[-2 maxFWx1/m/g],'Color',Colors(5,:),...
        'LineWidth',1,'LineStyle',Style(2));
else
    p1(5) = line([phimax1,phimax1],[-2 2],'Color',Colors(4,:),...
        'LineWidth',1,'LineStyle',Style(2));
end

if FrictionCase >0
    legend(p1(1:8),'F_{Bz}','F_{Wx}','\phi_s  \mu=0',...
     'max. Stabil. \mu=0','F_{Bz} \mu>0','F_{Wx} \mu>0',...
     '\phi_s \mu>0', 'max. Stabil. \mu>0',...
     'location','northeast','NumColumns',2);
else
    legend(p1(1:5),'F_{Bz}','F_{Wx}','\phi_s  \mu=0',...
          'max. Stabil. \mu=0', 'max. Stabil. \mu>0',...
          'location','northeast','NumColumns',3);
end    

grid on                                                                                                                 
ylabel('F in mg','FontSize',14)
xlabel('\phi in °','FontSize',14)
legend box off
h=title('Kräfte');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%% Funktionen

% Lagrangegleichung
function dY = dgl_Ladder01(t,Y,P1)
    % Y(1)- Position phi(t), 
    % Y(2)- Winkelgeschwindigkeit dphi(t)
    y  = Y(1);
    dy = Y(2);
    m     = P1.m;
    L     = P1.L;
    g     = P1.g;

    dY = [dy;...
          -3*g*cos(y)/2/L];
end

function dY = dgl_Ladder02(t,Y,P1)
    % Y(1)- Position phi(t), 
    % Y(2)- Winkelgeschwindigkeit dphi(t)
    y  = Y(1);
    dy = Y(2);
    m     = P1.m;
    L     = P1.L;
    g     = P1.g;
    mu    = P1.mu;
    dY = [dy;...
          -3/(2-mu^2)*(g/L*(1-mu^2).*cos(y)-mu.*sign(dy).*(dy.*dy-2*g.*sin(y)/L))];
end


function [value,isterminal,direction] = MyEvent(t,Y,P1)
    value = (0.0-rad2deg(Y(1))); 
    isterminal = 1;      % stop Integration
    direction = 0;       % negative Richtung
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


