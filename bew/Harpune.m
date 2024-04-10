% -------------------------------------------------------------------------
% Harpune.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Dynamik des Enterhakens auf Basis der
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

L       = 10.0;             % Länge Seil m
mS      = 10.0;             % Masse des Seils kg
mK      =  5.0;             % Masse des Körpers des Springers kg
mu      = mS/L;             % Massendichte in kg/m 
gamma   = mS/mK;

g       = 9.81;             % g in m/s 
z0      = 0;                % Anfangsgposition der Seils
dz0     = 10;               % Wurfgeschwindigkeit m/s
tmax    = dz0/g;            % Steigzeit freier Wurf
H       = dz0^2/2/g;        % Steighöhe freier Wurf
tspan   = linspace(0.0,tmax,100);

fprintf('\n ');
fprintf('\n L   = %8.2f m', L);
fprintf('\n mu  = %8.2f kg/m', mu);
fprintf('\n mS  = %8.2f kg/m', mS);
fprintf('\n mK  = %8.2f kg/m', mK);
fprintf('\n ');
fprintf('\n z0  = %8.2f m', z0);
fprintf('\n dz0 = %8.2f m/s', dz0);
fprintf('\n g   = %8.2f m/s^2', g);
fprintf('\n tmax= %8.2f m/s', tmax);
fprintf('\n ');


%% Berechnungen

% Anfangswerte
AB=[z0;dz0]; % AB für ode45

gamma  = [0.2, 1, 10.0];

gammastr(1,:) = "Freier Wurf";

for m= 1:3
    mK = gamma(m)*mS;
    gammastr(m+1,:) = strcat('\gamma = ',string(num2str(gamma(m),'% 5.2f')));
    % Numerische Lösung volle Lagrangegleichung
    opt=odeset('AbsTol',1.e-7,'RelTol',1.e-6);
    % Numerische Lösung LGL
    [t,Y]=ode45(@dgl_BungeeJumpRev,tspan,AB,opt,g,mK,mu); 
    for k=1:length(t)
        if Y(k,1) > L/2
           Y(k,1) = NaN;
           Y(k,2) = NaN;
        end
    end
    
    switch m
    case 1
        t1   = t;
        z1   = Y(:,1);
        dz1  = Y(:,2);
        ddz1 = g*(1+(mS.*z1.*(4*mK*L+2*mS*L-mS.*z1))./(2*(mS*L-mS.*z1+2*mK*L).^2));
   case 2
        t2   = t;
        z2   = Y(:,1);
        dz2  = Y(:,2);
        ddz2 = g*(1+(mS.*z2.*(4*mK*L+2*mS*L-mS.*z2))./(2*(mS*L-mS.*z2+2*mK*L).^2));
    case 3
        t3   = t;
        z3   = Y(:,1);
        dz3  = Y(:,2);
        ddz3 = g*(1+(mS.*z3.*(4*mK*L+2*mS*L-mS.*z3))./(2*(mS*L-mS.*z3+2*mK*L).^2));
    end
end

zf = dz0*t2-0.5*g*t2.^2;
dzf= dz0-g*t2;


%% Graphische Ausgabe

% Lösung Lagrangegleichung z(x)
fig=figure();
set(fig,'defaultAxesColorOrder',[Colors(3,:); Colors(2,:)]);
subplot(1,2,1)

h=title('Numerische Lösung Lagrange-Gl.');
hold on
lp(4) = plot(t2,abs(zf),'Color',Colors(4,:), 'LineWidth',2,'LineStyle',Style(1));
lp(1) = plot(t2,abs(z1),'Color',Colors(3,:), 'LineWidth',2,'LineStyle',Style(2));
lp(2) = plot(t2,abs(z2),'Color',Colors(3,:), 'LineWidth',2,'LineStyle',Style(1));
lp(3) = plot(t2,abs(z3),'Color',Colors(3,:), 'LineWidth',2,'LineStyle',Style(3));
axis([0 tmax 0 1.1*H]);
grid on
ylabel('Höhe z in m','FontSize',14)
xlabel('t in s','FontSize',14)
legend(lp, gammastr, 'location','south');
legend box off
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(1,2,2)
hold on;
rp(1) = plot(t1,abs(zf)-abs(zf),'Color',Colors(4,:), 'LineWidth',2,'LineStyle',Style(1));
rp(2) = plot(t1,abs(z1)-abs(zf),'Color',Colors(2,:), 'LineWidth',2,'LineStyle',Style(2));
rp(3) = plot(t1,abs(z2)-abs(zf),'Color',Colors(2,:), 'LineWidth',2,'LineStyle',Style(1));
rp(4) = plot(t1,abs(z3)-abs(zf),'Color',Colors(2,:), 'LineWidth',2,'LineStyle',Style(3));
ylabel('Abweichung v. freien Wurf in m','FontSize',14)
xlabel('t in s','FontSize',14)
axis([0 tmax -inf  0]);
legend(rp, gammastr, 'location','south');
legend box off
grid on
set(gca,'FontSize',16);

%%
% Berechnung Tauwurf

mS0 = 10;
mK = 7.5;
L  = 25;
mu = mS/L;
H = 5;
vf0 = sqrt(2*g*H);
v0  = linspace (vf0,2*vf0,20);

for k=1:5
 mS = k*mS0;
 mu(k) = mS/L;
 gammastr(k,:) = strcat('m_{Seil} = ',string(num2str(mS,'% 4.0f kg')));
 z_max(k,:) = mK/mu(k)*((1+3*mu(k)*v0.*v0/2/mK/g).^(1/3)-1);
end
figure();
hold on;

fprintf('\n Berechnung Tauwurf ');
fprintf('\n L   = %8.2f m', L);
fprintf('\n mu  = %8.2f kg/m - %4.2f kg/m ', mu(1), mu(4));
fprintf('\n mK  = %8.2f kg/m', mK);
fprintf('\n H   = %8.2f m', H);
fprintf('\n ');

for k=1:5
   h(k) = plot(v0,z_max(k,:),'Color',Colors(k+1,:), 'LineWidth',2);
end
line([0 3*vf0], [H H],'Color',Colors(9,:), 'LineWidth',2,'LineStyle',':');
axis([vf0 2*vf0 H/3 1.5*H]);
ylabel('Wurfhöhe in m','FontSize',14)
xlabel('Wurfgeschwindigkeit v_0 in m/s','FontSize',14)
legend(h, gammastr,'Location','Southeast');
legend box off;
grid on
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktioen

% Lagrange-Gleichung
function dY = dgl_BungeeJumpRev(t,Y,g,mK,mu)
    % Y(1)- Position z(t), 
    % Y(2)- Geschwindigkeit dz(t)
    dY    = [Y(2);...    
         -g - (mu./(mK+mu.*Y(1))).*Y(2).^2];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

