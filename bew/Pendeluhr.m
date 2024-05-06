% -------------------------------------------------------------------------
% Pendeluhr.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Pendeluhr
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", "-", "--", ":","-."];

%%
% Parameter: 

LR    = 0.150; % Laenge Regelschieber
RKS   = 0.055; % Radius Kreisscheibe
LP    = 0.400; % Laenge Aufhängung
g     = 9.81;

% Massen
mR    = 0.048;
mKS   = 0.156;
mP    = 0.016;
mKS   = 0.156;

%% 
kend = 4;
JP   = (1/3)*mP*LP*LP;
for k = 1:kend
  lr(k)   = 0.05 + (k-1)*0.01;
  for m = 1:95
    lks(m)   = 0.200 + (m-1)*0.002;
    % Traegheitsmomente bzgl. Aufhaengung, Regelschieber Parameter
    JKS(k,m)  = (1/2)*mKS*RKS*RKS + mKS*lks(m).*lks(m);
    JR(k,m)   = (1/12)*mR*LR*LR + mR*(lr(k) + LR/2).^2;
    Jges(k,m) = JP + JKS(k,m)+JR(k,m);
    % Schwerpunkt
    Mges = mR + mP +mKS;
    R(k,m)    = (mR*(lr(k)+LR/2)+mKS*lks(m)+mP*LP/2)/Mges;
    % Periodendauer
    limit = (lr(k)+LR);
    if lks(m)  > limit
        TP(k,m)  = 2*pi*sqrt(Jges(k,m)./(Mges*g*R(k,m)));
    else
        TP(k,m) = NaN;
    end
        
  end
end

figure()
hold on
for k=1:kend
    plot(lks*1000, TP(k,:),'Color',Colors(k+1,:),'LineWidth',2,...
        'LineStyle', Style(k));
    lgdstr(k,:)= strcat('\it l_R \rm= ',num2str(lr(k)*1000,'%4.0f mm'));
end
axis([256 258 0.998 1.002]);
legend(lgdstr,'location','southeast','FontSize',14,'FontName','times',...
      'NumColumns',2); 
legend box off
grid on
ylabel('\it T_P \rm in s','FontSize',14,'FontName','times');
xlabel('Scheibenposition \it l_{KS} \rm in mm','FontSize',14);
set(gca,'FontSize',14);

fig=figure();
ksel = 3;
kdtl = 9;
set(fig,'defaultAxesColorOrder',[Colors(ksel+1,:); Colors(kdtl,:)]);
hold on
yyaxis left
plot(lks*1000, TP(3,:),'Color',Colors(ksel+1,:),'LineWidth',2,...
        'LineStyle', Style(1));
lgdstr2(1,:)= lgdstr(1,:);
axis([240 340 0.85 1.15]);
ylabel('\it T_P \rm in s');

yyaxis right
plot(lks*1000, R(ksel,:)*1000,'Color',Colors(kdtl,:),'LineWidth',2,...
        'LineStyle', Style(2));
plot(lks*1000, Jges(ksel,:)*10000,'Color',Colors(kdtl,:),'LineWidth',2,...
        'LineStyle', Style(4));
% axis([200 350 0.0 1500]);
ylabel('\it R \rm  in mm, J_{ges} \rm in kg cm^2 ','FontSize',14, ...
       'FontName','times','color',Colors(kdtl,:));
hold off

legend('\it T_p','\it R(l_{KS})','\it J_{ges}(l_{KS})',...
      'location','east');
legend box off
grid on
xlabel('Scheibenposition \it l_{KS} \rm in mm','FontSize',14);
title(lgdstr(ksel,:));
set(gca,'FontSize',14);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

