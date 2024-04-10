% -------------------------------------------------------------------------
% Schaukel03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Schaukel mit Überschlag
% 
% Lösung der DGL für die Schaukel mit Kette und verhindertenm Überschlag
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Erwachsener
L       = 2.0;                      % Länge Schaukel
mg      = 100;                      % Gesamtmasse
g       = 9.81;                     % Erdbeschleunigung
tmax    = 4;
NPts    = 500;
t       = linspace(0,tmax,NPts);    % Zeitbereich
step    = tmax/NPts;
TSPAN   = 0:step:tmax;
kend    = 5;


figure(); % Trajektorien
hold on
axis([-L L -L L])
for k = 1:kend
    v0(k) = 6. + (k-1)*1;
    theta0      = 0;
    theta0dot   = v0(k)/L;  
    %Anfangswerte
    AB=[theta0;theta0dot];           % AB für ode45
    %Numerische Lösung
    Ereignis = @(T, Y) MyEventFunction(T, Y, g, L);
    options = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'Events',Ereignis);
    [T,Y, tE, YE, IE]=ode45(@(t,Y)dgl_schaukel(t, Y, g, L), ...
                                               TSPAN, AB,options); 
    x = L*sin(Y(:,1))-(k-1)*0.25;
    z = -L*cos(Y(:,1));
    sintheta = x/L;
    costheta = -z/L;
    zmax = v0(k)^2/2/g - L;         % Energiemaximum
    if k == 3 
       F_Tx = -L*Y(:,2).^2.*sin(Y(:,1)) - g.*sin(Y(:,1)).*cos(Y(:,1));
       theta_Tx = Y(:,1);
    end
    h(k) = plot(x,z,'Color', Colors(k,:), 'LineWidth',2);
    iend = length(x);
    AB2=[L*sin(Y(iend,1));-L*Y(iend,2)*sin(Y(iend,1));...
         -L*cos(Y(iend,1));L*Y(iend,2)*cos(Y(iend,1));];     % AB für ode45
    Ereignis2 = @(T, Y) MyEventFunction2(T, Y, L);
    options2 = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'Events',Ereignis2);
    [T,Y]=ode45(@(t,Y)dgl_fall(t, Y, g),TSPAN, AB2,options2); 
    % Hier noch Abfangen für Überschlag bei L^2 = x^2+z^2, 
    % Frage ob ich durch "Aufschaukeln",
    % so viel Energie gewinnen kann, dass Zentrifugalkräfte ausreichen für 
    % Überschlag von Kettenschaukel
    xf = Y(:,1)-(k-1)*0.25;
    zf = Y(:,3);
    if v0(k)< sqrt(5*g*L) 
       plot(xf,zf,'Color', Colors(k,:), 'LineWidth',2);
    else
       plot(-(x+(k-1)*2*0.25),z,'Color', Colors(k,:), 'LineWidth',2);
    end
    x1 = [-L L];
    z1 = [zmax zmax];
    p1 = line(x1 ,z1);
    p1.Color = Colors(k,:);
    p1.LineWidth = 1;
    x2 = [0-(k-1)*0.25 x(end)];
    z2 = [0 z(end)];
    p2 = line(x2 ,z2);
    p2.Color = Colors(k,:);
    p2.LineWidth = 1;
    p2.LineStyle = ':';
    x3 = [0-(k-1)*0.25 0-(k-1)*0.25];
    z3 = [0 -L];
    p3 = line(x3 ,z3);
    p3.Color = Colors(k,:);
    p3.LineWidth = 1;
    p3.LineStyle = ':';
    p4 = line([-L L],[0 0]);
    p4.Color = Colors(15,:);  p4.LineWidth = 2;
    p4.LineStyle = '-';
    strtemp      = num2str(v0(k),4);   
    lgdstr(k,:)  = string(sprintf(' v_0 = %s m/s ', strtemp));
end
axis([-L L -L L])
grid on
% legend(h,lgdstr,'location','bestoutside');
xlabel('$x$','interpreter','latex','FontSize',16)
ylabel('$z$','interpreter','latex','FontSize',16)
h=title('Schaukel');
set(h,'FontSize',12,'FontWeight','normal'); 
axis square;
set(gca,'FontSize',16);


figure(); % Zwangskräfte
hold on
axis([-L L -L L])
for k = 1:kend
    v0(k) = 6. + (k-1)*1;
    theta0      = 0;
    theta0dot   = v0(k)/L;  
    % Anfangswerte
    AB=[theta0;theta0dot];       % AB für ode45
    % Numerische Lösung
    Ereignis = @(T, Y) MyEventFunction(T, Y, g, L);
    options = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'Events',Ereignis);
    [T,Y, tE, YE, IE]=ode45(@(t,Y)dgl_schaukel(t, Y, g, L), ...
                                               TSPAN, AB,options); 
    x = L*sin(Y(:,1))-(k-1)*0.25;
    z = -L*cos(Y(:,1));
    sintheta = x/L;
    costheta = -z/L;
    zmax = v0(k)^2/2/g - L; %Energiemaximum
    F_Tx = -L*Y(:,2).^2.*sin(Y(:,1)) - g.*sin(Y(:,1)).*cos(Y(:,1));
    theta_Tx = Y(:,1);
    h(k)=plot(rad2deg(theta_Tx),F_Tx,'Color',Colors(k,:), 'LineWidth',2);
end
axis([0 180 min(F_Tx)*1.1 g])
legend(h,lgdstr,'location','southeast');
legend box off;
xlabel('$\theta$','interpreter','latex','FontSize',16)
ylabel('$F_{T,x}$','interpreter','latex','FontSize',16)
grid on
h=title('Kettenspannung');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
function dY = dgl_schaukel(t, Y, g,L)
    % Y(1)-Winkel
    %Y(2)-Winkelgeschwindigkeit
    dY = [Y(2);...
         -g*sin(Y(1))/L];
end

function dY = dgl_fall(t, Y, g )
    % Y(1)-Winkel
    % Y(2)-Winkelgeschwindigkeit
    dY = [Y(2);0;...
          Y(4);-g];
end


% Ereignisfunktion
function [value,isterminal,direction] = MyEventFunction(t, Y, g, L)
    % Locate the time when tension passes through 
    % zero 
    % and stop integration.  
    value = -L*Y(2)^2*sin(Y(1))-g*cos(Y(1))*sin(Y(1));  % detect force = 0
    isterminal = 1;   % stop the integration
    direction  = 0;   % all direction
end

  
% Ereignisfunktion
function [value,isterminal,direction] = MyEventFunction2(t, Y, L)
    % Locate the time when fall passes through 
    % pendelum length 
    % and stop integration.  
    value = L^2-Y(1)^2-Y(3)^2;  % detect pendelum length 
    isterminal = 1;   % stop the integration
    direction  = 0;   % all direction
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

