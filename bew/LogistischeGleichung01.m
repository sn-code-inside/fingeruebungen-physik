% -------------------------------------------------------------------------
% LogistischeGleichung01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Getriebener Nichtlinearer Oszillator/Pendel
% Logistische Gleichung
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% eine Trajektorie
figure();

NPts = 75;
Gamma1 = 2.8;       
% Gamma1 = 3.5; %alternativ              
x(1)    = 0.2;
x(2)    = Gamma1*x(1)*(1-x(1));
line([x(1) x(1)], [0 x(2)],'Color',Colors(3,:), 'LineWidth',1);
hold on
line([x(1) x(2)], [x(2) x(2)],'Color',Colors(3,:), 'LineWidth',1);
for k = 2:NPts+1
  x(k+1) = Gamma1*x(k)*(1-x(k));    
  line([x(k) x(k)], [x(k) x(k+1)]);
  line([x(k) x(k+1)], [x(k+1) x(k+1)]);
end
xcoord = 0:0.01:1;
yccord = Gamma1*xcoord.*(1-xcoord);
plot(xcoord,yccord, 'color',Colors(4,:),...
           'LineStyle', Style(3),'Linewidth',1);  
plot(xcoord,xcoord, 'color',Colors(4,:),...
           'LineStyle', Style(4),'Linewidth',1);  
grid on
axis square;
axis equal;
axis([0 1 0 1])
xlabel('$x_k$','interpreter','latex','FontSize',14)
ylabel('$x_{k+1}$','interpreter','latex','FontSize',14)
Gammastr = cat(2,'\Gamma = ', num2str(Gamma1,5));
ht = text(0.1,0.9, Gammastr); 
set(ht,'FontSize',12);
grid on

% verschiedene Trajektorien
NPts   = 50;
kt     = 1:NPts+1;
NPlots = 5;
xt = zeros(NPlots,NPts);
for m = 1:NPlots
    Gamma(m) = 0.75*m;
    xt(m,1)    = x(1);
    xt(m,2)    = Gamma(m)*xt(m,1)*(1-xt(m,1));
    GammaStr(m) = string(cat(2,'$\Gamma = $', num2str(Gamma(m),5)));
end
figure();
hold on
for m = 1:NPlots
    for k = 2:NPts
      xt(m,k+1) = Gamma(m)*xt(m,k)*(1-xt(m,k));    
    end
    hp(m) = plot(kt(:),xt(m,:),'color',Colors(7-m,:),'LineWidth',1,...
            'LineStyle',Style(1));
end
xlabel('Zyklus \it{k}','FontSize',14)
ylabel('Population \it{x}','FontSize',14)
axis([1 NPts 0 1]);
grid on;
hlp = legend(hp, GammaStr ,'interpreter','latex'); 
set(hlp,'FontSize',12);
set(hlp,'location', 'bestoutside');
legend box off;
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

