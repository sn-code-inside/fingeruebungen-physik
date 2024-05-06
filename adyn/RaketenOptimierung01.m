%%-----------------------------------------------------------------------%%
%
% Raketenoptimierung.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
%
%--------------------------------------------------------------------------
%
% Programm optimiert Parameter einer Mehrstufenrakete nach der
% Extremwertmethode
% Als Beispiel wird die Saturn V mit Ihren Maximalwerten genommen, einmal
% als zweistufige Rakte und einmal als dreistufige.
% 
%--------------------------------------------------------------------------

%%

clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%%
% Daten Saturn V 3-stufig 

sigmaS  = [0.05	0.06 0.06];                 % Strukturkoeffizienten
vGS     = [2761	4561 4549];                 % Gasgeschwindiglkeiten in m/s
m0S     = [2828500 458700+217250 217250];   % Untermassen in kg
mLS     = 100000;                           % Nutzlast in kg
muLS    = mLS/m0S(1);

 % Optimierungsbereich Nutzlastverhältnis
muL      = linspace(0.0,0.06,120);         
 % Vorbereitung
for i=2:3
  Q(i)  = (vGS(i)-vGS(i-1))/(2*vGS(i-1)*sigmaS(i));
  P(i)  = vGS(i)*sigmaS(i-1)/vGS(i-1)/sigmaS(i);
end
Q2 = Q(2);
Q3 = Q(3);
P2 = P(2);
P3 = P(3);
 % Lösung der kubischen Gleichung
for k=1:length(muL)
    b = 2*muL(k)*(2*Q3+Q2*P3);
    c = 4*Q3*muL(k).*muL(k)*(Q3+Q2*P3);
    d = P2*P3*P3*muL(k).*muL(k);
    p = [1 -b c -d];
    r(:,k) = roots(p);
end
mu3 = r(1,:);
mu2 = Q2*mu3+sqrt(Q2^2*mu3.^2+P2*mu3);

% Optimiertes Antriebsvermögen für gegebenes µ_L
Deltav=vGS(1)*log(1./(sigmaS(1)+mu2))+ ...
       vGS(2)*log(1./(sigmaS(2)+mu3./mu2)) +...
       vGS(3)*log(1./(sigmaS(3)+muLS./mu3));

% Alternative Berechung für µ_3 über Numerisch-Graphische Lösung
% muL = 0.004;  %Vorgegebenes µ_L
% mu3 = linspace(0.02,0.06,400); %Suchbereich µ_3
% W2  = sqrt(Q2^2.*mu3.^2+P2.*mu3);
% F   =  Q3.*muL +...
%         sqrt(Q3^2.*muL.^2+P3.*muL.*Q2.*mu3 +P3*muL.*W2);
% figure(2)
% plot(mu3,F,)
% hold on
% plot(mu3,mu3)
% ylabel('µ_3, F(µ_3)','FontSize',14)
% xlabel('µ_L','FontSize',14)
% set(gca,'FontSize',16);



%%
% Daten Saturn V 2-stufig angenommen
sigmaS  = [0.05	0.06 0.06];           % Strukturkoeffizienten
vGS     = [2761	4561 4549];           % Gasgeschwindiglkeiten in m/s
m0S     = [2828500-117250	675950-117250	0];  % Untermassen in kg
mLS     = 100000;                                % Nutzlast in kg

 % Vorbereitung
for i=2:3
  Q(i)  = (vGS(i)-vGS(i-1))/(2*vGS(i-1)*sigmaS(i));
  P(i)  = vGS(i)*sigmaS(i-1)/vGS(i-1)/sigmaS(i);
end
muLS = mLS/m0S(1);
mu22 = Q(2).*muL+sqrt((Q(2).*muL).^2+P(2).*muL);

% Optimiertes Antriebsvermögen für gegebenes µ_L (analytische Formel)
Deltav2= vGS(1)*log(1./(sigmaS(1)+mu22))+ ...
         vGS(2)*log(1./(sigmaS(2)+muLS./mu22));

   

%% 
% Graphische Ausgabe

newcolors = [Colors(3,:)
             Colors(2,:)
             Colors(4,:)
             Colors(8,:)];       
figure('Name','Zwei-Drei-Stufen Optimierung')
set(groot,'defaultAxesColorOrder',newcolors)
subplot(1,2,1)
yyaxis right
% ax = uiaxes;
h(1)=plot (muL, mu2,'color',Colors(2,:),'Linewidth',2);
hold on
h(2)=plot (muL, mu22,'color',Colors(2,:),'Linewidth',2,...
                     'LineStyle', Style(2));
ylabel('µ_2','FontSize',14)
yyaxis left
h(3)=plot (muL, mu3,'color',Colors(3,:),'Linewidth',2);
hold on
x = [muLS muLS];
y = [min(mu3) max(mu3)];
h(4)=line(x,y,'color',Colors(4,:),'Linewidth',2,'LineStyle', Style(3));
grid on;
legend(h,'3 Stufen \mu_2','2 Stufen \mu_2','3 Stufen \mu_3',...
         'Saturn V real','location','northwest','numColumns',1)
legend box off
ylabel('µ_3','FontSize',14)
xlabel('µ_L','FontSize',14)
set(gca,'FontSize',16);
set(gca,'defaultAxesColorOrder','remove')


subplot(1,2,2)
% plot (muL, Deltav/1000,'color',Colors(2,:),'Linewidth',2);
semilogy(muL, Deltav/1000,'color',Colors(2,:),'Linewidth',2);
hold on
% plot (muL, Deltav2/1000,'color',Colors(2,:),'Linewidth',2,...
%      'LineStyle', Style(3));
semilogy (muL, Deltav2/1000,'color',Colors(2,:),'Linewidth',2,...
     'LineStyle', Style(3));
x = [muLS muLS];
yminAx = 5;
ymaxAx = 15;
y = [yminAx ymaxAx];
pl = line(x,y,'color',Colors(4,:),'Linewidth',2,'LineStyle', Style(3));
grid on;
ylabel('Antriebsvermögen km/s','FontSize',14)
ylim([yminAx ymaxAx]);
xlabel('µ_L','FontSize',14)
legend('3 Stufen ','2 Stufen','Saturn V real',...
       'location','southwest','NumColumns',1)
legend box off
set(gca,'FontSize',16);

