% -------------------------------------------------------------------------
% LambertProblem02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Parameter der Ellipsenbahn 
% für das Lambert-Problem für einen Marsflug durch ein Bisektionsverfahren.
%
% Benutzt LambertSolver2.m.
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
muS  = 2.95479E-04;      % [AE^3/d^2] für Sonne

%% Beginn Rechnung
r1    = 1.0000;          % in AE
r2    = 1.525;           % in AE
gamma = deg2rad(75);
chord = sqrt(r1^2+r2^2-2*r1*r2*cos(gamma));
s     = (chord+r1+r2)/2;
TOF   = 115;   % in d
mu    = muS;


%% Graphische Lösung------------------------------------------------------% 
%
amin  = s/2;
alpham  = 2*asin(sqrt(s/2/amin));
alphamd = rad2deg(alpham);
betam   = 2*asin(sqrt((s-chord)/2/amin));
betamd  = rad2deg(betam);
tm      = sqrt(amin^3/mu)*(alpham-betam-(sin(alpham)-sin(betam)));
a       = linspace(amin,5*amin,1000); 
alpha1  = 2*asin(sqrt(s/2./a));
alpha2  = 2*pi-2*asin(sqrt(s/2./a));
for k=1:length(alpha1)
   if ~isreal(alpha1(k))
       alpha1(k) = NaN;
   end
   if ~isreal(alpha2(k))
       alpha2(k) = NaN;
   end
end
beta   = 2*asin(sqrt((s-chord)/2./a));
tpar   = sqrt(2)/3 * sqrt(s^3/mu) * (1-sqrt((s-chord)^3/s^3));
tF1    = sqrt(a.^3/mu).*(alpha1-beta-(sin(alpha1)-sin(beta)));
tF2    = sqrt(a.^3/mu).*(alpha2-beta-(sin(alpha2)-sin(beta)));
% Graphik
figure(1)
h(1) = plot(a, tF1,'color',Colors(2,:),'Linewidth', 2);
hold on
h(2) = plot(a, tF2,'color',Colors(3,:),'Linewidth', 2);
h(3) = line([0 1.23],[365 365],'color',Colors(3,:),'Linewidth',1,... 
'LineStyle',Style(3));
h(4)= line([0 1.23],[115 115],'color',Colors(2,:),'Linewidth',1,... 
'LineStyle',Style(3));
h(5)= line([0 2.0],[tpar tpar],'color',Colors(9,:),'Linewidth',1,... 
'LineStyle',Style(3));
h(6)= line([0.0 amin],[tm tm],'color',Colors(4,:),'Linewidth',1,... 
'LineStyle',Style(3));
line([1.23 1.23],[0 365],'color',Colors(3,:),'Linewidth',1,... 
'LineStyle',Style(3));
h(7)= line([1.23 1.23],[0 115],'color',Colors(2,:),'Linewidth',1,... 
'LineStyle',Style(3));
h(8)= line([amin amin],[0 tm],'color',Colors(4,:),'Linewidth',1,... 
'LineStyle',Style(3));
grid on
ylabel('\tau(a) in Tagen')
xlabel('a in AE')
axis([0 2 0 1000]); 
lgd = legend(h,"\tau (\alpha)","\tau (\alpha ')","t_F (oberer Ast)",...
             "t_F (oberer Ast)",...
             "t_{Par}","t_m","a (gesucht)","a_{min}",...
             'location','bestoutside');
legend box off;
set(lgd,'FontSize',12,'FontWeight','normal');
ttl=title('Lambert Problem \tau(a)');
set(ttl,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',16);

%% Iteration -------------------------------------------------------------% 

[aiter,amin,amax,alpha_n,beta_n,tau_n,tm,tpar,nend] = ...
          LambertSolver2(r1, r2, gamma, TOF, mu);

aend = aiter(nend);
% Graphik
figure(2)
hp(1) = plot(a,tF1,'color',Colors(2,:),'Linewidth',2,'LineStyle',Style(2));
hold on
hp(2) = plot(aiter(1:nend-1),tau_n(2:nend),'o-','LineWidth',1);
for k=1:nend-1
  if k <10
      text(aiter(k)-0.001,tau_n(k+1)-0.5,num2str(k,2),'Color',Colors(4,:));
  else
      text(aiter(k)-0.001,tau_n(k+1)+0.5,num2str(k,2),'Color',Colors(4,:)); 
  end
end
plot(aiter(nend),tau_n(nend+1),'d-','LineWidth',2, 'MarkerSize',6,...
    'color', Colors(2,:),'MarkerFaceColor',Colors(2,:))
hp(3)= line([0 aiter(nend)],[TOF TOF],'color',Colors(2,:),'Linewidth',1,... 
     'LineStyle',Style(3));
hp(4)= line([aiter(nend) aiter(nend)],[0 TOF],'color',Colors(2,:),...
     'Linewidth',1,... 
'LineStyle',Style(3));
grid on
ylabel('\tau(a) in Tagen')
xlabel('a in AE')
axis([0.9*aend 1.1*aend 0.9*TOF 1.1*TOF]); 
lgd = legend(hp,"\tau (\alpha)","Bisektionsiteration","t_F", ...
             "a (gesucht)",'location','northeast');
legend box off;
set(lgd,'FontSize',12,'FontWeight','normal');
ttl=title('Lambert Problem - Bisektionsverfahren');
set(ttl,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',16);

% PrintOut
ttlprint = "Iteration Lambert-Problem für Mars-Flug";
fprintf('\n %s', ttlprint);
fprintf('\n')
fprintf('\n    n   |   a_n  (AE)   |   alpha_n   |   beta_n     |  tau(a_n)(d) | tau(a_n)-TOF (d)  |   e_n   \n\n');
for k=1:nend
psi = (alpha_n(k)-beta_n(k))/2;
phi = (alpha_n(k)+beta_n(k))/2;
e_n(k) = sqrt((r2-r1)^2/(2*aiter(k)*sin(psi))^2+(cos(phi))^2);

fprintf('  %3u   |   %7.4f     |   %6.3f    |   %6.3f     |  %6.1f      |   %10.3f      |  %7.4f \n',...
        k, aiter(k), alpha_n(k), beta_n(k), tau_n(k+1),  tau_n(k+1)-TOF, e_n(k));
end
fprintf('\n');

%-------------------------------------------------------------------------%
% Ende Programm
% ------------------------------------------------------------------------%

% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
