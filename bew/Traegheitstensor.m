% -------------------------------------------------------------------------
% Traegheitstensor.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnung des Traegheitstensors 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%
% ellipso.m calculates an ellipsoid's inertia tensor, & mass
% numerically, also plots it in 3d
% Ellipsoid: (x-x0)^2/a^2+(y-y0)^2/b^2+(z-z)^2/c^2=1
% uses Simpson's rule for 3d integration
a=3;b=2;c=1;                                   % semimajor axes
rho=1/8;                                       % density 
ax=-a;                                         % x lower limit
bx=a; Nx=35; dx=(bx-ax)/(Nx-1);                % x upper limit, points, spacing
x=[ax:dx:bx];                                  % x grid
Ny=35;Nz=35;                                   % y,z points
ie=7;                                          % matrix elements calculated + volume
for k=2:Nx                                     % x loop - begin at 2 so ay ~= 0
    ay=-b*sqrt(1-(x(k)/a)^2);                  % y lower limit,
    by=b*sqrt(1-(x(k)/a)^2);                   % y upper limit,
    if real(by) ~=0                            % check if array exists
      dy=(by-ay)/(Ny-1);                       % y spacing
      y=[ay:dy:by];                            % y grid
      for j=2:Ny                               % y loop - begin at 2 so az ~= 0
         az=-c*sqrt(1-(x(k)/a)^2-(y(j)/b)^2);  % z lower limit
         bz=c*sqrt(1-(x(k)/a)^2-(y(j)/b)^2);   % z upper limit
            if real(bz) ~= 0                   % check z array exists
              dz=(bz-az)/(Nz-1);               % z spacing
              z=[az:dz:bz];                    % z grid
              for i=1:Nz                       % z loop
                for m=1:ie
                  fz(m,i)=rho*inert_el2(m,x(k),y(j),z(i));
                end
              end                              % end z loop
            end                                % end 2nd if
         for m=1:ie
            fy(m,j)=simp(fz(m,:),dz);          % Simpson rule
         end
      end                                      % end y loop
    end                                        % end 1st if
    for m=1:ie
        fx(m,k)=simp(fy(m,:),dy);              % Simpson rule
    end
end                                            % end x loop
% finally integrate over the x coord to get the moments
for m=1:ie
    ff=simp(fx(m,:),dx);                       % Simpson rule
       if     m==1 Ixx=ff;
       elseif m==2 Ixy=ff;
       elseif m==3 Ixz=ff;
       elseif m==4 Iyy=ff;
       elseif m==5 Iyz=ff;
       elseif m==6 Izz=ff;
       elseif m==7 Mass=ff;
       end
% fprintf('m= %2i, The integral is %4.3f\n',m,ff)  
end
Iyx=Ixy; Izx=Ixz; Izy=Iyz;                    % use symmetry for the rest
A=[[Ixx,Ixy,Ixz];[Iyx,Iyy,Iyz];[Izx,Izy,Izz]];% inertia tensor
disp(['Semimajor axes (m): a,b,c=',rat(a),' ',rat(b),' ',rat(c)])
disp 'Inertia Tensor in kgm^2'
disp(rats(A))                        % display A in string fraction form
M=rho*4*pi*a*b*c/3;                  % actual ellipsoid mass
p_e=(Mass-M)*100/M;                  % error on the mass
str1=cat(2,'density=',num2str(rho,3),'kg/m^3, Mass=',num2str(Mass,3),...
           ' kg, % mass error=',num2str(p_e,3));
disp(str1)
% draw the ellipsoid centered at x0,y0,z0, & semimajor axes a,b,c
x0=0; y0=0; z0=0;N=50;
[x,y,z]=ellipsoid(x0,y0,z0,a,b,c,N); % uses 50 mesh points
h=mesh(x,y,z,'EdgeColor',[0.0 0.0 0.9]);
axis equal, grid on, box on, view ([-0.1 0.5 0.2])
xlabel('x','FontSize',14),ylabel('y','FontSize',14),zlabel('z','FontSize',14)
title(['Ellipsoid: a=',rat(a),', b=',rat(b),', c=',rat(c),', \rho=',...
      num2str(rho,3),', M=',num2str(M,3),', I_{xx}=',num2str(A(1,1),2),...
      ', I_{yy}=',num2str(A(2,2),2),', I_{zz}=',num2str(A(3,3),2)],...
      'FontSize',13)
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


% inert_el2.m
function inerel2=inert_el2(m,x,y,z)    
  % Inertia function integrands in cartesian coords
  % if m=7 it does the volume integrand
  if     m==1 inerel2=y^2+z^2;
  elseif m==2 inerel2=-x*y;
  elseif m==3 inerel2=-x*z;
  elseif m==4 inerel2=x^2+z^2;
  elseif m==5 inerel2=-y*z;
  elseif m==6 inerel2=x^2+y^2;
  elseif m==7 inerel2=1.0;
  else
    disp ' only 7 integrands are needed '
    return
  end
end

% simp.m
function simpu=simp(f,inc)   
  % Simpson's rule for numerical integration
  % f is an odd array of evaluated functions in steps inc
  ip=length(f);       %must be an odd number
  s1=sum(f(2:2:ip-1));%sums all even terms
  s2=sum(f(3:2:ip-2));%sums all odd term does not include f(1) and f(ip)
  simpu=(4.*s1+2.*s2+f(1)+f(ip))*inc/3.0;%finally add f(1) and f(ip)
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
