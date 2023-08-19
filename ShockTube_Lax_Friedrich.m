%% 1D Shock Tube using Lax-Friedrich and Adjustable Time Stepping
%close all
clearvars -except pressure_plot_save u_plot_save T_plot_save e_plot_save rho_plot_save
clc
%% Initialization of Parameters
N       = 2000;             % Number of grid points
%gamma   = 1.4 ;
%gamma = 1.313; %methane
Cv = 1709.;
Cp = 2232.;
Rgas = Cp - Cv;
gamma = Cp / Cv;



dia = 1.15; %pipe inner diameter
rough = 0.05 * 1e-3;
lambda = 1 / (-2*log10(rough/(3.7 * dia)))^2;
%lambda = 0.1; %fanning friction coefficient

%lambda = 1;
%
%lambda = 0.0;

area = dia.^2/4*pi;

save_ = 1;
load_ = 0;

shocktube = 1;

endTime = 3600 * 24 * 3;             % Physical End Time
%endTime = 60/433;
%endTime = 400/433;

dtplot = endTime / 100;
%dt_pplot = 1;
dt_pplot = endTime / 10;
CFL     = 0.3;
%% Grid Initialization
%x       = linspace(0,1,N+1) ;
%x       = linspace(0, 433 * endTime * 0.1,N) ;
x = linspace(0, 230000, N);
%x = linspace(0, 1000000, N);
xc = x;
%x = linspace(0,100, N + 1);
dx = xc(2) - xc(1);
%xc      = 0.5 * ( x(1:N) + x(2:N+1) ) ;
%xc(1)   = 0 ;               % Xmin
%xc(N)   = 1 ;               % Xmax
%xc(N)   = 100 ;               % Xmax
time    = 0 ;
%% Initial Conidtions
%denistyR    = 0.125 ;       densityL    = 1.0 ;
%pressureR   = 0.1   ;       pressureL   = 1.0 ;

%

if (shocktube)

    pressureR   = 0.1 * 100000   ;       pressureL   = 165 * 100000 ;
    %%"nordstream" conditions

    %pressureR   = 165 * 100000   ;       pressureL   = 400 * 100000 ;
    
    T0 = 280;
    
    %densityR    = 0.1*0.72 ;       densityL    = 0.72 * 164 ;
    densityR    = pressureR / (T0 * Rgas) ;       densityL    = pressureL / (T0 * Rgas);
    
    cl = sqrt(gamma*pressureL/densityL);
    umid = 2 / (1 + gamma) * cl;
    
    fact = 1. - 0.5 * (gamma-1) * umid / cl;
    
    rhomid = densityL * fact ^ (2. / (gamma - 1));
    
    pmid = pressureL * fact ^ (2. * gamma / (gamma-1));
    
    rho     = zeros(N,1) ;
    p       = zeros(N,1) ;
    u       = zeros(N,1) ;
    
    
    %PV=nRT => PV ~ rho*T
    
    for i =  1:N
        %if i<=N/2
        if (i < N)
            rho(i)  = densityL  ;
            p(i)    = pressureL ;
        else
            rho(i)  = densityR  ;
            p(i)    = pressureR ;
        end
    end
    
    
    
    if (load_)
        load('workspace.mat');
        time = 0;
        inc  = ceil(200 / dx);
        rho(inc:N/2) = rho(1:N/2 - inc + 1);
        rho(1:inc - 1) = densityL;
        rho(N/2+1:end) = densityR;
    
        p(inc:N/2) = p(1:N/2 - inc + 1);
        p(1:inc - 1) = pressureL;
        p(N/2+1:end) = pressureR;
    
        u(inc:N/2) = u(1:N/2 - inc + 1);
        u(1:inc - 1) = 0;
        u(N/2+1:end) = 0;
    
        %u(401:end) = 0.;
        %u(1:end) = 0.;
    end

else

    T0 = 280;
    pressureL   = 165 * 100000;
    densityL    = pressureL / (T0 * Rgas);
    %pressureS = 330 * 100000;
    %pressureL / densityL^gamma == pressureS / densityS^gamma

    densityS = 2 * densityL;
    pressureS = pressureL * (densityS / densityL) ^ gamma;

    pressureR = pressureL;
    densityR = densityL;
    
    cl = sqrt(gamma*pressureL/densityL);

    %(2 * ncell + 1) * dx = 1 => (1/dx - 1) / 2

    %ncell = ceil( (1. / dx - 1) / 2);
    ncell = 100;

    rho = ones(N, 1) * densityL;
    rho(N / 2 - ncell:N/2 + ncell) = densityS;

    p = ones(N, 1) * pressureL;
    p(N / 2 - ncell:N/2 + ncell) = pressureS;

    u = zeros(N, 1);


end

u0 = u;
p0 = p;
rho0 = rho;

T = p ./ (rho .* Rgas);

% for i =  1:N
%     if i<=N - 1
%     %if i<=N
%         rho(i)  = densityL  ;
%         p(i)    = pressureL ;
%     else
%         rho(i)  = denistyR  ;
%         p(i)    = pressureR ;
%     end
% end


e   = p/(gamma-1) + 0.5*rho.*u.*u ;

%p = (gamma - 1) * (e  - 0.5 * rho * u * u)


%%
new_rho = rho ;
new_u   = u   ;
new_e   = e   ;
new_T = T;

totmass = [];
mass_right = [];
mass_left = [];
e_right = [];
e_left = [];
plot_times = [];
tot_e = [];

pressure_plot = [];
e_plot = [];
u_plot = [];
rho_plot = [];
T_plot = [];

plot_time = time;

pplot_time = time;

e(N-1) = e(N-2) * 0.5;

while time <= endTime
%     
%     for i=2:N-1
%         p(i)    = (gamma-1)*(e(i) - 0.5*rho(i)*u(i)*u(i)) ;
%     end

    %p(2:N-1) = (gamma-1)*(e(2:N-1) - 0.5*rho(2:N-1).*u(2:N-1).*u(2:N-1)) ;
    %T(2:N-1) = p(2:N-1) ./ (rho(2:N-1) * Rgas);
    p = (gamma-1)*(e - 0.5*rho.*u.*u) ;
    T = p ./ (rho * Rgas);

    %p_  / (new_rho(i) * Cv * T(i) * (gamma - 1)

    a       = sqrt(gamma*p./rho) ;
    lamda   = max(a) ;
    max_vel = max(u) ;
    
    dt      = CFL/(max_vel+lamda) * dx ;  % adjustable Time Step
    time    = time + dt ;

    fprintf("time: %f\n", time);
    
    for i=2:N-1
        

        mom_R       = rho(i+1)*u(i+1) ;     rho_R = rho(i+1) ;      u_R = u(i+1) ;      p_R = p(i+1) ;
        mom_P       = rho(i)*u(i)     ; 	rho_P = rho(i)   ;      u_P = u(i)   ;      p_P = p(i)   ;
        mom_L       = rho(i-1)*u(i-1) ;    	rho_L = rho(i-1) ;      u_L = u(i-1) ;      p_L = p(i-1) ;

        lambda_P = abs(u(i)) + sqrt(gamma * p(i) / rho(i));
        lambda_L = abs(u(i - 1)) + sqrt(gamma * p(i - 1) / rho(i - 1));
        lambda_R = abs(u(i + 1)) + sqrt(gamma * p(i + 1) / rho(i + 1));
        lambda_L = max(lambda_L, lambda_P);
        lambda_R = max(lambda_R, lambda_P);
        
        vel_flux_R  = rho_R*u_R*u_R +p_R ;    e_R = e(i+1) ;
        vel_flux_P  = rho_P*u_P*u_P +p_P ;    e_P = e(i)   ;
        vel_flux_L  = rho_L*u_L*u_L +p_L ;    e_L = e(i-1) ;
        
        energy_flux_R   = u_R * ( e_R + p_R ) ;
        energy_flux_P   = u_P * ( e_P + p_P ) ;
        energy_flux_L   = u_L * ( e_L + p_L ) ;
        
        rho_fluxR   = 0.5*( mom_P + mom_R ) -0.5*lambda_R*( rho_R - rho_P ) ;
        rho_fluxL   = 0.5*( mom_P + mom_L ) -0.5*lambda_L*( rho_P - rho_L ) ;
        
        vel_fluxR   = 0.5*( vel_flux_P + vel_flux_R ) -0.5*lambda_R*( mom_R - mom_P ) ;
        vel_fluxL   = 0.5*( vel_flux_P + vel_flux_L ) -0.5*lambda_L*( mom_P - mom_L ) ;
        
        energy_fluxR    = 0.5*( energy_flux_P + energy_flux_R ) -0.5*lambda_R*( e_R - e_P );
        energy_fluxL    = 0.5*( energy_flux_P + energy_flux_L ) -0.5*lambda_L*( e_P - e_L );

        if (i == N - 1)
        
        %if (0)
            energy_fluxR = energy_flux_P;
            if (0)
                vel_fluxR = vel_flux_P;
            end
            rho_fluxR = mom_P;
        end
        



        vel_flux    = mom_P - dt/dx * ( vel_fluxR - vel_fluxL ) ;
        
        new_rho(i)  = rho_P - dt/dx * ( rho_fluxR - rho_fluxL ) ;
        
        %vel_flux = mom_P - dt/dx * ( vel_fluxR - vel_fluxL ) - dt*fric + dt*0;
        
        new_u(i)    = vel_flux/new_rho(i) ;
        new_e(i)    = e_P - dt/dx * ( energy_fluxR - energy_fluxL ) ;

        if (i == 2)
            new_rho(1) = rho(1) - dt/dx * rho_fluxL;
            mom1 = rho(1)*u(1) - dt/dx * (vel_fluxL - (0*rho(1)*u(1).^2 +p(1)));
            new_u(1) = mom1 / new_rho(1);
            new_e(1) = e(1) - dt/dx * (energy_fluxL - 0*(u(1) * ( e(1) + p(1) )));
        end
% 
%         if (i == N - 1)
%             new_rho(N) = rho(N) + dt/dx * rho_fluxR;
%             new_u(N) = u(N) + dt/dx * vel_fluxR;
%             new_e(N) = e(N) + dt/dx * vel_fluxR;
%         end

        
        %if (shocktube)
        if (0)
            if (i <= N/2) 
                u_fric = (2 * u_P) / (1+sqrt(1+(16*dt*u_P * lambda) / dia));
        %fric = 0;
            else
                u_fric = u_P;
            end
        else
            %u_fric = (-dia + sqrt(dia*(dia + 16 * dt * u_P* lambda))) / (8 * dt * lambda);
            %fric = rho_P * sign(u_P)*u_P^2 / 8 * lambda / dia * 4.; %dimension Pa/m
            u_fric = (2 * u_P) / (1+sqrt(1+(16*dt*u_P * lambda) / dia));
        end

        %dU_fric = -dt * fric / new_rho(i);
        dU_fric = u_fric - u_P;
        dE_fric = new_rho(i) * ((new_u(i) + dU_fric)^2 - new_u(i)^2) / 2; %kinetic energy lost to friction

        new_u(i) = new_u(i) + dU_fric;
        %new_e(i) = new_e(i) + dE_fric;

        
    end



%     vel_flux_1 = p(1);
% 
%     new_u(1) = new_u(1) + dt/dx * vel_flux_1;
% 
%     rho_flux_N = rho(N) * u(N);
%     vel_flux_N = rho(N) * u(N).^2 + p(N);
%     energy_flux_N = u(N) * (e(N) + p(N));
% 
%     new_rho(N) = new_rho(N) - dt/dx * rho_flux_N;
%     new_u(N) = new_u(N) - dt/dx * vel_flux_N;
%     new_e(N) = new_e(N) - dt/dx * energy_flux_N;



    %new_u(1) = new_u(2);

    %new_u(N) = new_u(N - 1);
    %new_rho(N) = new_rho(N - 1);
    %new_e(N) = new_e(N - 1);

    %fprintf("totmass: %f\n", sum(new_rho));
    if (time >= plot_time)
        plot_time = plot_time + dtplot;
        plot_times = [plot_times time];
        totmass = [totmass sum(new_rho)];
        mass_right = [mass_right sum(new_rho(N/2:end))];
        mass_left = [mass_left sum(new_rho(1:N/2-1))];
        tot_e = [tot_e sum(new_e)];
        e_left = [e_left sum(new_e(1:N/2-1))];
        e_right = [e_right sum(new_e(N/2:end))];
    end

    if (time >= pplot_time) 
        pplot_time = pplot_time + dt_pplot;
        pressure_plot = [pressure_plot p];
        e_plot = [e_plot e];
        u_plot = [u_plot u];
        rho_plot = [rho_plot rho];
        T_plot = [T_plot T];
    end

    %new_u(1) = new_u(2);
    %new_rho(1) = new_rho(2);
    %new_e(1) = new_e(2);

    new_u(N) = new_u(N - 1);
    %new_u(N) = 0;
    new_rho(N) = new_rho(N - 1);
    %new_rho(N) = densityR;
    new_e(N) = pressureR / (gamma-1) + 0.5*new_rho(N - 1) .* new_u(N - 1) .* new_u(N - 1);
    %new_e(N) = new_e(N - 1);
    
    %new_Tn = p(N - 1) / rho(N - 1) / Rgas;
    %new_rho(N) = pressureR / (new_Tn * Rgas);
    %T(2:N-1) = p(2:N-1) ./ (rho(2:N-1) * Rgas);

    %new_e(N/2+1:end) = pressureR / (gamma-1) + 0.5*new_rho(N/2+1:end) .* new_u(N/2+1:end) .* new_u(N/2+1:end);
    %new_u(N/2) = sqrt(gamma * p(N/2) / rho(N/2));
    
     rho     = new_rho ;
     u       = new_u ;
     e       = new_e ;

     %rho(1:400)     = new_rho(1:400) ;
     %u(1:400)       = new_u(1:400) ;
     %e(1:400)       = new_e(1:400) ;

     %rho(N) = rho(N - 1);
     %rho(N) = densityR;
     %u(N) = u(N - 1);
     %e(N) = u(N)^2 * rho(N)/2 + 1 / (gamma - 1) * pressureR;
     
     
end

if (save_)
    save("workspace.mat");
    rho_plot_save = rho_plot;
    e_plot_save = e_plot;
    T_plot_save = T_plot;
    pressure_plot_save = pressure_plot;
    u_plot_save = u_plot;
end

% pressure = dlmread('pressure.dat') ;
% density  = dlmread('density.dat')  ;
% velocity = dlmread('velocity.dat') ;
% 
% figure(1)
% hold on
% %plot(density(:,1),density(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
% plot(xc, rho,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
% xlabel(' x ','FontSize',18,'FontWeight','bold');
% ylabel(' Density ','FontSize',18,'FontWeight','bold');
% legend('Lax Friedrich','Location','northeast','FontSize',16,'FontWeight','bold');
% %print(gcf,'Density.jpg','-dpng','-r300');
% 
% figure(2)
% hold on
% %plot(pressure(:,1),pressure(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
% plot(xc, p,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
% xlabel(' x ','FontSize',18,'FontWeight','bold');
% ylabel(' Pressure ','FontSize',18,'FontWeight','bold');
% legend('Lax Friedrich','Location','northeast','FontSize',16,'FontWeight','bold');
% %print(gcf,'Pressure.jpg','-dpng','-r300');
% 
% figure(3)
% hold on
% %plot(velocity(:,1),velocity(:,2),'r--','MarkerSize',12,'MarkerIndices',1:10:length(velocity),'LineWidth',2);
% plot(xc, u,'k','MarkerSize',12,'MarkerIndices',1:10:141,'LineWidth',2);
% xlabel(' x ','FontSize',18,'FontWeight','bold');
% ylabel(' Velocity ','FontSize',18,'FontWeight','bold');
% legend('Lax Friedrich','Location','south','FontSize',16,'FontWeight','bold');
%print(gcf,'Velocity.jpg','-dpng','-r300');