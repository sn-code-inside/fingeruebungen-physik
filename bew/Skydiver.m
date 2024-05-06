% -------------------------------------------------------------------------
% Skydiver.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Freier Fall aus großer Höhe (Baumgartner-Weltrekord).
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];
Marker = ['o','d','o','s','+'];

% Variablen: 
m      = 130;       % in kg
kappa  = 1.40e-4;   % in 1/m
g      = 9.81;      % in m/s^2
A      = 1.0;       % in 1/m^2 Annahme Kopf nach unten
cw     = 0.2;
rho0   = 1.225;     % in 1/m^3
eta0   = 0.5*A*cw*rho0/m;

% Höhe
z0   = 38969;       % in m
z    = linspace(7500,z0,1000);

% Integrationszeit
tmax = 300;         % in s
tspan=[0.0,tmax];

h0 = 38969;
q  = 2*eta0*exp(-kappa*z)/kappa;
q0 = 2*eta0/kappa*exp(-kappa*z0);
C  = - expint(q0);
v  = sqrt(-2*g/kappa.*exp(-q).*(expint(q)+C));

% Messwerte
MW(1,:)  = [33.446, 310];
MW(2,:)  = [27.883, 377];
MW(3,:)  = [22.961, 290];
MW(4,:)  = [07.619, 079];
MW(5,:)  = [02.567, 053];

ISATable = ReadExcelTable(42,'F');
hz   = ISATable(:,1)/1000;
rhoz = ISATable(:,2);
cwz  = ISATable(:,5);
Az   = ISATable(:,6);

z2    = linspace(0,40000,1000);
rho1 = rho0*exp(-kappa*z2);
figure()
plot(z2/1000,rho1,'Linewidth',1,'Color',Colors(2,:),'LineStyle',Style{1});
semilogy(z2/1000,rho1,'Linewidth',1,'Color',Colors(2,:),'LineStyle',Style{1});
hold on;
semilogy(hz,rhoz,'Linewidth',1,'Color',Colors(4,:),'LineStyle',Style{1});
axis([0 40  0.001 1.25]);
legend('Barometr. Höhenformel','ISA Tabelle 1',...
'FontSize',14, 'location', 'best');
set(gca,'XDir','reverse');
grid on
legend box off
xlabel('Höhe \it z \rm in km ','FontSize',14);
ylabel('Luftdichte in kg/m^3 ','FontSize',14);
set(gca,'FontSize',16)   

% Numerische Berechnung
h0 = 39000;
m  = 100;
[tout,zout] = falling_guy(h0, hz, g, rhoz, Az, cwz,m);

% Graphische Darstellung
figure()
plot(z/1000,v,'Linewidth',2,'Color',Colors(3,:),'LineStyle',Style{1});
axis([0 40 0 400]);
hold on
plot(zout(:,1)/1000, abs(zout(:,2)),'Linewidth',2,'Color',Colors(4,:));
axis([0 40 0 400]);
set(gca,'XDir','reverse');
grid on
xlabel('Höhe \it z \rm in km ','FontSize',14);
ylabel('Geschwindigkeit \it v \rm in m/s ','FontSize',14);
for k = 1:5 
    plot(MW(k,1),MW(k,2),Marker(k),'MarkerSize',8,'MarkerEdgeColor',...
            Colors(k+1,:),'MarkerFaceColor',Colors(k+1,:),'LineWidth',2);
end
legend('Analyt. Näherung','Numer. Berechnung', ...
       'Mach 1','Maximale v', 'Mach 1', 'Vorb. Schirm',...
       'Öffnung Schirm','FontSize',14, 'location', 'best');
legend box off;   
grid on;
set(gca,'FontSize',16)   
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%% Funktion zur numerischen Berechnung

function [tout,zout] = falling_guy(h0, hz, g, rhoz, Az, cwz,m);
tstart = 0;
tfinal = 500;
AB=[h0;0];       % AB für ode45
% Mit Plot
% options = odeset('Events',@events,'OutputFcn',@odeplot,'OutputSel',2);
% figure();
% set(gca,'xlim',[0 300],'ylim',[0 40000]);
% set(gca,'xlim',[0 300],'ylim',[-500 0]);
% hold on;

% Ohne Plot
options = odeset('AbsTol',1.e-9,'RelTol',1.e-5,'Events',@events);

tout = tstart;
teout = [];
zout = AB';

zeout = [];
veout = [];
ieout = [];
for k = 1:36
    % Solve until the first terminal event.
    hi = 40-k;
    eta = 0.5*Az(hi)*cwz(hi)*rhoz(hi)/m;
    P1 = hz(hi)*1000;
    [t,y,te,ye,ie] = ode23(@(t,y)[y(2);-g+eta*y(2)*y(2)],...
                     [tstart tfinal],AB,options);
    if ~ishold
        hold on
    end
    % Accumulate output.  This could be passed out as output arguments.
    nt = length(t);
    tout = [tout; t(2:nt)];
    zout = [zout; y(2:nt,:)];
    vend = y(end,2);
    ieout = [ieout; ie];
    tend(k)=t(end);
    % Set the new initial conditions
    y0(1) = hz(hi)*1000;
    y0(2) = vend;
    AB = [y0(1),y0(2)];
    tstart = t(nt);
end


% Nested function
    function [value,isterminal,direction] = events(t,y)
        % Locate the time when height passes through zero in a decreasing direction
        % and stop integration.  Here we use a nested function to avoid
        % passing the additional parameter P1 as an input argument.
        value = y(1) - P1;     % detect height = 0
        isterminal = 1;   % stop the integration
        direction = 0;   % negative direction
    end
end 

    

% %% Funktionen
% %-------------------------------------------------------------------------
% % DGL Freier Fall aus großer Höhe
% % In  : rhoz. czwP
% % Out : dzdt
% function dY = dgl_fall(t, Y,eta)
% % Calculate parameters
% g    = 9.81;
% eta  = 0.02;
% %  Es  muss  ein  Spaltenvektor  zurückgegeben  werden 
% dY     =  zeros(2,1); 
% dY(1)  =  Y(2);
% dY(2)  = -g+eta*Y(2)*Y(2);
% end



%% Einlesen der ISA Standardtabelle
% In  : ZE Zeilenend (Zahl)  SE SpaltenEnde (String)
% Out : ISATable     ..  Werte  (Tabelle)
%-------------------------------------------------------------------------
 
%_________________________________________________________________________
%
function ISATable = ReadExcelTable(ZE,SE)
    fname = 'ISA01.xlsx';
    fields = strcat('A1:',SE,num2str(ZE));
    fid=fopen(fname,'r');
    if fid ==1 
        disp('File open not successful');
    else
        ISATable = xlsread(fname,1,fields);
    end
    closeresult =fclose(fid);
    if closeresult ==0
    %     disp('Data successfully loaded');
    else
        disp('Dataload not successful');
    end
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


