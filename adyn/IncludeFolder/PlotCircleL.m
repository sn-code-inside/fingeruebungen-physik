function PlotCircleL (xM,yM,rC,Col,LW, LS)
    u=linspace(0,360,360);
    nx=zeros(360);
    ny=zeros(360);
    Nmax = 1; 
    for k=1:Nmax
        nx= rC*k*cosd(u)/Nmax+xM;
        ny= rC*k*sind(u)/Nmax+yM;
        plot(nx,ny,'LineWidth',LW,'Color',Col,'LineStyle',LS);
    end
end

