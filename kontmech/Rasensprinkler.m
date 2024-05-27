% -------------------------------------------------------------------------
% Rasensprinkler.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Charakteristika eines rotierenden Rasensprinklers
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

dD  = 0.004; AD= pi/4*dD^2;      % Austrittsquerschnitt in m^2 
din = 0.0125; Ain = pi/4*din^2;  % Eintrittsquerschnitt in m^2 
R   = 0.30;               % Halbe Armlänge in m
Q   = 0.1e-3;             % Fluss in m^3/s  (1/10 l pro s)
phi = deg2rad(10);        % Winkel zur Tangente
rho = 1e3;                % Dichte in kg/m^3
omeg= 0.5;                % Umdrehungszahl in 1/s = 30 min^-1
pL  = 1013e2;             % Luftdruck in Pa
eta = Q*R*rho*(Q*cos(phi)/2/AD - R);  % in Nms 
p   = pL + 0.5*rho*(Q^2/AD^2-omeg^2*R^2-Q^2/Ain^2);  % in Pa 

% omega als Funktion von Q und AD


figure(1)
eta = 10;
Qvec   = linspace(0.01e-3,0.5e-3,40);
dDvec  = linspace(0.002,0.006,5);
ADvec  = pi/4*dDvec.^2;

for k=1:length(ADvec)
    omega(k,:) = Qvec*cos(phi)./ADvec(k)/2/R^2./(1+eta./Qvec/rho/R^2);
end
omegvec= omega;
Deltap = omega;
hold on
axis([0 0.5 0 60]);
for k=1:length(ADvec)
    plot(Qvec(:)*1000,omega(k,:)*60,'color',Colors(k+1,:),'linewidth',2);
    lgdstr(k,:) = strcat(string(num2str(dDvec(k)*1000,4)),' mm');
end
xlabel('Fluss Q in l/s ','FontSize',14)
ylabel('Umdrehung in 1/min','FontSize',14)
grid on
set(gca,'FontSize',14);
streta =strcat(num2str(eta,3),' Nms )');
str = "Rotation als Funktion Durchfluss Q und Düsendurchmesser d_D  (Reibung \eta = ";
str = strcat(str,streta);
h2 = title(str);
legend(lgdstr, 'location','northwest');
legend box off
set(h2,'FontWeight','normal','FontSize',12);
h2 = title(str);

figure(2)
eta = 10;
dD     = 0.005;
AD     = pi/4*dD^2;
etavec= linspace(5,25,5)/60;
for k=1:length(etavec)
    omegvec(k,:) = Qvec*cos(phi)./AD/2/R^2./(1+etavec(k)./Qvec/rho/R^2);
    Deltap(k,:)  = 0.5*rho*(Qvec.^2/AD^2-Qvec.^2/Ain^2-omegvec(k,:).^2*R^2);
    Deltap(k,:)  = 1 + Deltap(k,:)/pL;
end
hold on
% axis([0 1 0 3]);
for k=1:length(etavec)
    plot(Qvec(:)*1000,Deltap(k,:),'color',Colors(k+1,:),'linewidth',2);
    lgdstr(k,:) = strcat(string(num2str(etavec(k)*60,2)),' Nms');
end
xlabel('Fluss Q in l/s ','FontSize',14)
ylabel('Überdruck p/p_L','FontSize',14)
grid on
set(gca,'FontSize',14);
strD  = strcat(num2str(dD*1000,2),' mm');
strin = strcat(num2str(dD*1000,2),' mm');
str  = strcat("Druck als Funktion Durchfluss Q und Reibungszahl \eta ( d_D = ",strD);
str  = strcat(str,', d_{in} =  ');
str  = strcat(str, strin);
str  = strcat(str, ' )');
h2=title(str);
legend(lgdstr, 'location','northwest');
legend box off
set(h2,'FontWeight','normal','FontSize',12);
h2 = title(str);

