% -------------------------------------------------------------------------
% FallenderStab01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Kräfte beim fallenden Stab 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Numerische Lösung

%Parameter

L       =  0.2;             % Länge Stab in  m
H       =  L/2;             % Schwerpunktlage bei homogenen Stab  m
m       =  0.01;            % Masse des Stabs in kg
mu1     =  0.10;            % Reibungskoeffizient 1 
mu2     =  0.30;            % Reibungskoeffizient 2 
mu3     =  0.60;            % Reibungskoeffizient 3 
JP      =  m*L^2/3;         % Trägheitsmoment Fusspunkt   
JS      =  m*L^2/12;        % Trägheitsmoment Schwerpunkt
g       =  9.81;            % g in m/s^2 
theta0  =  deg2rad(01.0);   % Anfangswinkel rad
dtheta0 =  0;               % Anfangsgeschwindigkeit rad/s
tmax    =  sqrt(2*H/g);     % Fallzeit freier Fall
tspan   =  linspace(0.0,1,90);

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

%%
% Kräfteberechnung als Funktion des Winkels 
theta = linspace(deg2rad(1),deg2rad(90),100);
omega02=m*g*H/JP;
FH   = m*H*omega02*sin(theta).*(3*cos(theta)-2*cos(theta0));
FN   = m*g-m*H*omega02*(1-3*cos(theta).*cos(theta)+2*cos(theta).*cos(theta0));
% Bestimmung Slip-Winkel thetaS
myfun = @(x,c1,c2) (sin(x).*(9*cos(x)-6*cos(c2))./...
                   (1+9*cos(x).*cos(x)-6*cos(x).*cos(c2))) - c1;  
                 % parameterized function               
thetaC = acosd(2*cos(theta0)/3);                
for k=1:3
    fun = @(x) myfun(x,mu(k),theta0);    % function of x alone
    thetaS(k) = rad2deg(fzero(fun,theta0));
    if (thetaS(k) < 0) || (thetaS(k) > thetaC)
        thetaS(k) = thetaC;
    end
end

fig=figure();
set(fig,'defaultAxesColorOrder',[Colors(3,:); Colors(2,:)]);
hold on
for k=1:3 
    fp(k) = plot(rad2deg(theta),1000*mu(k)*FN,'Color',Colors(k,:));
end
fp(4) = plot(rad2deg(theta),1000*FH,'Color',Colors(4,:));
set(fp,'LineWidth',2,'LineStyle',Style(1));
for k=1:3 
fp(k+4) = line([thetaS(k) thetaS(k)],[-1000 1000],...
        'Color',Colors(k,:),'LineStyle',Style(2));
end
fp(8) = line([thetaC thetaC],[-1000 1000],...
        'Color',Colors(8,:),'LineStyle',Style(3));
set(fp,'LineWidth',2);
maxY = 1100*max(mu)*max(FN);
axis([0 90 -maxY maxY]);  
grid on
xlabel('Winkel \theta  in °','FontSize',14)
ylabel('Kraft in mN','FontSize',14)
legend(fp,strcat('\muF_N @ \mu_1=',num2str(mu1,3)),...
          strcat('\muF_N @ \mu_2=',num2str(mu2,3)),...
          strcat('\muF_N @ \mu_3=',num2str(mu3,3)),'F_H',...
       '\theta_S @ \mu_1', '\theta_S @ \mu_2','\theta_S @ \mu_3',...
       '\theta_C ','location','eastoutside');
legend box off
ht=title('Kräfte beim fallenden Stab als Funktion des Winkels bei festen Fußpunkt');
set(ht,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


