% -------------------------------------------------------------------------
% JoJo02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet Dynamik des akademischen Jo-Jos mit Handbewegung
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

L       =  1.000;           % Länge Stab in  m
Ri      =  0.009;           % Innenradius in m
Ra      =  0.045;           % Aussenradius in m
g       =  9.81;            % g in m/s^2 
m       =  0.040;           % Masse in kg 
x0      =  0;               % Anfangsposition
dx0     =  0.0;             % Anfangsgeschwindigkeit in m/s
tmax1   =  5;               % max Zeitspanne für Berechnung


fprintf('\n ');
fprintf('\n L   = %8.2f m',  L);
fprintf('\n m   = %8.2f kg', m);
fprintf('\n Ri  = %8.2f mm', Ri*1000);
fprintf('\n Ra  = %8.2f mm', Ra*1000);
fprintf('\n g   = %8.2f m/s^2', g);
fprintf('\n ');


%% Berechnungen

% Anfangswerte
AB1=[x0,dx0]; % AB für ode45

%%
% Parameter

P1.g      = g;
P1.Ra     = Ra;
P1.Ri     = Ri;
P1.L      = L;
P1.tstart = 0;

% Umlenkung
omegaB   =  2*sqrt(g*L)/Ra;
tmax2    =  2*pi/omegaB;
tspan1   =  linspace(0.0,tmax1,200);

% Optionen für ODE45
options1 = odeset('AbsTol',1.e-7,'RelTol',1.e-5,'events',@MyEvent1);
options2 = odeset('AbsTol',1.e-7,'RelTol',1.e-5,'events',@MyEvent2);
options3 = odeset('AbsTol',1.e-7,'RelTol',1.e-5,'events',@MyEvent3);

% Felder
z     = NaN(4,200);
z2    = NaN(4,200);
z3    = NaN(4,200);
zr    = NaN(4,200);
phi   = NaN(4,200);
v2    = NaN(4,200);
td    = NaN(4,200);
td2   = NaN(4,200);
td3   = NaN(4,200);
FT2   = NaN(4,200);

% Parameter Set für Schleife
aHup    = [0,0,g/2, 2*g/3];
aHdown  = [0,0,g/3, 1*g/2];
tHup    = [0,0,0.2,0.15];
tHdown  = [0,0,0.2,0.1];
eta     = [0,0.25,0.25,0.25];
kend    = length(aHup);

% Schleife für verschieden Parameter
for k = 1:kend
P1.aHup   = -aHup(k);           % Stärke Impuls up
P1.aHdown = +aHdown(k);         % Stärke Impuls down
P1.tHup   = tHup(k);
P1.tHdown = tHdown(k);
P1.eta    = eta(k);             % Dämpfung

Para_Str1(k,:) = strcat(strcat(' |aH| =  ', ...
                       num2str(aHup(k)/g,'% 5.2f')),' g ; ');
Para_Str2(k,:) = strcat(Para_Str1(k,:),strcat( ' \eta =  ', ...
                       num2str(P1.eta,'% 5.2f')),' 1/s ');
% Berechnung Abrollphase
[t, Y, TE, YE, IE] = ode45(@dgl_JoJo1,tspan1, AB1, options1,P1); 
    for k1=1:length(t) 
        td(k,k1) = t(k1); 
        z(k,k1)  = Y(k1,1);
        vend(k) = YE(2);
        tfall(k)= TE; 
    end
    for k1=1:length(z(k,:))
        zr(k,k1) = z(k,length(z(k,:))+1-k1);      
    end
% Berechnung Umkehrphase
    AB2=[+pi/2,-vend(k)/Ri]; % AB für ode45
    options2 = odeset('AbsTol',1.e-7,'RelTol',1.e-5,'events',@MyEvent2);
    tspan2   =  linspace(tfall(k),tfall(k)+tmax2,200);
    [t, Y, TE, YE, IE] = ode45(@dgl_JoJo2,tspan2, AB2,options2,P1); 
    for k1=1:length(t) 
        td2(k,k1) = t(k1); 
        phi(k,k1) = Y(k1,1);
        z2(k,k1)  = cos(Y(k1,1))*Ri;
        tturn(k)  = TE; 
        v2(k,k1)  = Ri*Y(k1,2)*sin(Y(k1,1));
    end
% Berechnung Aufrollphase
    AB3=[z2(k,length(t)),v2(k,length(t))]; % AB für ode45
    tspan3   =  linspace(tturn(k),tmax1+tturn(k),200);
    P1.tstart = tspan3(1);
    [t, Y, TE, YE, IE] = ode45(@dgl_JoJo3,tspan3, AB3, options3,P1); 
    for k1=1:length(t) 
        td3(k,k1)  = t(k1); 
        z3(k,k1)   = Y(k1,1);
    end

end
%% 
% Graphische Ausgabe

figure()
% Fallhöhe über t
subplot(1,2,1)
hold on
for k=1:4  
    p(k) = plot(td(k,:), z(k,:), 'Color',Colors(k,:),...
            'LineWidth',2,'LineStyle',Style(1));
    plot(td2(k,:), L+z2(k,:),'Color',Colors(k,:),...
            'LineWidth',2,'LineStyle',Style(1));
    plot(td3(k,:), L-z3(k,:),'Color',Colors(k,:),...
            'LineWidth',2,'LineStyle',Style(1));
end
ylabel('Fallhöhe h in m','FontSize',14)
axis([0, 2.5*tfall(1), -0.2, 1.2]);
axis ij
xlabel('t in s','FontSize',14)
grid on
h=title('Bewegtes Jo-Jo');
set(h,'FontSize',14,'FontWeight','normal'); 
legend(p, Para_Str2,'location','southeast','numcolumns',1);
legend box off
set(gca,'FontSize',16);

% Beschleunigung der Hand in Vielefachen von g
subplot(1,2,2)
hold on
for k=1:4   
    rp(k) = plot(td(k,:), aHdown(k)*(1-heaviside(td(k,:)-tHdown(k)))/g,...
            'Color',Colors(k,:), 'LineWidth',2,'LineStyle',Style(1));
     x = [0 0]; y = [0 aHdown(k)/g];
     line(x,y,'Color',Colors(k,:), 'LineWidth',2,'LineStyle',Style(1))
     plot(td3(k,:), -aHup(k)*(1-heaviside(td3(k,:)-tturn(k)-tHup(k)))/g,...
            'Color',Colors(k,:), 'LineWidth',2,'LineStyle', Style(1));
     x = [tturn(k) tturn(k)]; y = [0 -aHup(k)/g];
     line(x,y,'Color',Colors(k,:), 'LineWidth',2,'LineStyle',Style(1))
end
ylabel('a_H in g','FontSize',14)
axis([0, 2.5*tfall(1), -1, 1]);
axis ij
xlabel('t in s','FontSize',14)
grid on
h=title('Jo-Jo Beschleuingung Hand');
set(h,'FontSize',14,'FontWeight','normal'); 
legend(rp, Para_Str2,'location','southeast','numcolumns',1);
legend box off
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Functions

% Jo-Jo in Fallphase
function dY = dgl_JoJo1(t,Y,P1)
% Y(1)-  z(t), 
% Y(2)-  Geschwindigkeit dz(t)
rx    = P1.Ri; % akademisches Jo-Jo
dY    = [Y(2);...    
         (P1.g + P1.aHdown*(1-heaviside(t-P1.tstart-P1.tHdown)))/...
             (1+P1.Ra^2/2/P1.Ri^2) - P1.eta*Y(2)];
end


% Jo-Jo in Umkehrphase
function dY = dgl_JoJo2(t,Y,P1)
% Y(1)- Winkel phi(t), 
% Y(2)- Winkelgeschwindigkeit omega(t)
dY    = [Y(2);...    
         -2*sin(Y(1))*P1.g*P1.Ri/P1.Ra/P1.Ra];
end

 
% Jo-Jo in Aufwärtsphasephase
function dY = dgl_JoJo3(t,Y,P1)
% Y(1)-  z(t), 
% Y(2)-  Geschwindigkeit dz(t)
rx    = P1.Ri;
dY    = [Y(2);...    
         (-P1.g - P1.aHup*(1-heaviside(t-P1.tstart-P1.tHup)))/...
            (1+P1.Ra^2/2/P1.Ri^2) - P1.eta*Y(2)];
end



% Ereignisfunktionen

function [value,isterminal,direction] = MyEvent1(t,Y,P1)
% Ereignisfunktion bis vollständigen Abrollen
    value = (P1.L-Y(1)); % detect length of line
    isterminal = 1;      % stop the integration
    direction = 0;       % negative direction
end

function [value,isterminal,direction] = MyEvent2(t,Y,P1)
%Ereignisfunktion bis Eintreten der vollständigen Unkehr
    value = (-pi/2-Y(1)); % detect angle pi/2
    isterminal = 1;       % stop the integration
    direction = 0;        % negative direction
end

function [value,isterminal,direction] = MyEvent3(t,Y,P1)
% Ereignisfunktion bis vollständigen Aufrollen oder Geschwindigkeit 0
% thetaC
    value1 = Y(2);          % detect velocity 0
    value2 = P1.L-Y(1) ;    % detect position 0
    value  = value1*value2;
    isterminal = 1;      % stop the integration
    direction = 0;       % negative direction
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


