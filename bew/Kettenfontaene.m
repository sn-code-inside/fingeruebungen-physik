% -------------------------------------------------------------------------
% Kettenfontaene.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Z-Verlauf der Kettenfontäne
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
s = linspace(0,5,1000);
alpha   =   0.12;
beta    =   0.11;

H       =   [1.0,1.5,2.0,2.5,3];  
gamma = alpha*H/(1-alpha-beta);
for k=1:length(gamma)
    z(:,k) = (1-beta)/alpha * gamma(k) - abs(gamma(k)-s(:));
end

figure()
for k= 1:length(gamma)
     p(k) = plot(s,z(:,k),'Color',Colors(k,:),...
            'LineWidth',2);
     hold on
     parastr(k,:)  = string(strcat('H = ',num2str(H(k),'%4.1f m')));
end
H       =   2.5;                 %  
alpha   =   [0.12,0.05,0.1,0.15,0.2];
beta    = 0.11;
gamma = alpha*H./(1-alpha-beta);
for k=1:length(gamma)
    zalpha(:,k) = (1-beta)/alpha(k)*gamma(k) - abs(gamma(k)-s(:));
end
for k= 1:length(gamma)
     p(k+length(gamma)) = plot(-s,zalpha(:,k),'Color',Colors(4,:),...
            'LineWidth',2,'LineStyle', Style(k));
     hold on
     parastr(length(gamma)+k,:)  = string(strcat('alpha = ',num2str(alpha(k),'%4.2f')));
end
grid on
ymin = 0;
ymax = round(1.1*max(z(:,end)));
axis ([-s(end)+1,s(end)-1,ymin,ymax]);
ylabel('Höhe z in m','FontSize',14)
xlabel('Kettenlänge s in m','FontSize',14)
h=title('Kettenfontäne ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend([p],parastr, 'location','eastoutside',...
        'NumColumns',1);
legend box off
set(gca,'FontSize',16);

alpha = 0.22;
beta = 0.11;
H= linspace(0,4,100);
alpha   =   [0.12,0.05,0.1,0.15,0.2];
beta    =   [0.10,0.10,0.2,0.2,0.0];
for k=1:length(alpha)
  zmax(:,k) = (1-beta(k))*H(:)/(1-alpha(k)-beta(k))-H(:);
end

figure()
for k= 1:length(gamma)
     p2(k) = plot(H,zmax(:,k),'Color',Colors(k,:),...
            'LineWidth',1);
     hold on
     parastr2a(k,:)  = string(strcat('alpha = ',num2str(alpha(k),'%4.2f - ')));
     parastr2b(k,:)  = string(strcat('beta  = ',num2str(beta(k), '%4.2f ')));
     parastr2(k,:)  = strcat(parastr2a(k,:) ,parastr2b(k,:));
end
grid on
ymin = 0;
ymax = round(1.1*max(zmax(:,end)));
axis ([0,H(end),ymin,ymax]);
ylabel('Höhe \Delta z in m','FontSize',14)
xlabel('Höhe H in m','FontSize',14)
h=title('Kettenfontäne ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p2,parastr2, 'location','eastoutside',...
        'NumColumns',1);
legend box off
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


