% -------------------------------------------------------------------------
% MeteoritenEinschlag.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Einschlag eines Meteoriten aus großer Höhe.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":", "-."];
Marker = ['o','d','o','s','+'];


%% Numerischen Berechnung

% Integrationszeit
tmax = 20;         % in s
tsteps = 1000;
% Konstanten
RE     = 6370000;   % Erdradius in m
kappa  = 1.30e-4;   % in 1/m
g      = 9.81;      % in m/s^2
cw     = 0.45;
rho0   = 1.225;     % Dichte Atmosphäre am Bodenin kg/m^3
% Höhe 
z0   = 120000;      % in m
kend  = 6;
RM    = [0.01, 0.01, 0.01, 0.001, 0.005, 0.025];
for k = 1:kend 
    % Anfangsgeschwindigkeit
    if k < 4 
        v0 = k*15000; % in m/s  
    else
        v0 = 2*15000;         % in m/s  
    end  
    R = RM(k);
    C = {' v0 = ',num2str(v0/1000), ' km/s; ','R = ',num2str((R*1000),'%02d'),' mm'};
    lgdstr(k,:) = string(strjoin(C));
    % AB für ode45
    AB=[z0;-v0];    
    tspan=linspace(0.0,tmax,tsteps);
    rhoM   = 7900;      % Dichte Meteorid (Eisen) in kg/m^3
    A      = pi*R^2;    % in 1/m^2 
    m      = 4*pi*rhoM*R^3/3;
    eta    = 0.5*rho0*cw*A/m;
    options = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'Events',@Ereignis);
    [tout,zout] = ode45(@(t,Y)dgl_fall(t, Y, eta, kappa),tspan,AB,options);
    h(k,:) = zout(:,1)/1000;
    v(k,:) = abs(zout(:,2)/1000);
    temp = diff(zout(:,2)/g);
    temp = [temp;0];
    a(k,:) = temp(:);
 end

%% Graphische Darstellung

figure()
subplot (1,2,1)
hold on
for k = 1:kend
       if k < 4 
        kstyle = 1;
        kcolor = k+1;
    else
        kstyle = k;
        kcolor = 3;     
    end  
    plot(h(k,:),v(k,:),'Linewidth',2,...
      'Color',Colors(kcolor,:),'LineStyle',Style{kstyle});
end
axis([0 120 0 50]);
set(gca,'XDir','reverse');
grid on
xlabel('Höhe \it z \rm in km ','FontSize',14);
ylabel('Geschwindigkeit \it v \rm in km/s ','FontSize',14);
% legend(lgdstr,'location','best');
% legend box off;   
grid on;
set(gca,'FontSize',16)   


subplot (1,2,2)
hold on
for k = 1:kend
    if k < 4 
        kstyle = 1;
        kcolor = k+1;
    else
        kstyle = k;
        kcolor = 3;     
    end  
    plot(h(k,:),a(k,:),'Linewidth',2,...
      'Color',Colors(kcolor,:),'LineStyle',Style{kstyle});
end
axis([0 120 0 100]);
set(gca,'XDir','reverse');
grid on
xlabel('Höhe \it z \rm in km ','FontSize',14);
ylabel('Beschleunigung \it a \rm  in \it g ','FontSize',14);
legend(lgdstr,'location','best');
legend box off;   
grid on;
set(gca,'FontSize',16)   
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen
%-------------------------------------------------------------------------
% Ereignisfunktion
function [value,isterminal,direction] = Ereignis(t,y)
    value = y(1);     % detect height = 0
    isterminal = 1;   % stop the integration
    direction = 0;    % negative direction
end
  
%-------------------------------------------------------------------------
% DGL Freier Fall aus großer Höhe
function dY = dgl_fall(t, Y, eta, kappa)
    % Calculate parameters
    g    = 9.81;
    RE   = 6370000;
    dY     =  zeros(2,1); 
    q      = (RE/(RE+Y(1)))^2;
    dY(1)  =  Y(2);
    dY(2)  = -g*q+eta*exp(-Y(1)*kappa)*Y(2)*Y(2);
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


