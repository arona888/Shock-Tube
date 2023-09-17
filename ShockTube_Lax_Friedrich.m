%% 1D Shock Tube using Lax-Friedrich and Adjustable Time Stepping
%close all
clearvars -except pressure_plot_save u_plot_save T_plot_save e_plot_save rho_plot_save
clc
%% Initialization of Parameters

Cv = 1709.; %heat capacity of methane at constant volume (J/(K*Kg))
Cp = 2232.; %heat capacity of methane at constant pressure (J/(K*Kg))
Rgas = Cp - Cv;
gamma = Cp / Cv;

T0 = 285; %ambient temperature (K)
p0_ambient = 8 * 100000; %ambient pressure (Pa)
p0_pipe = 165 * 100000; %ns1 initial pipe pressure (Pa)
%p0_pipe = 103 * 100000; %ns2 initial pipe pressure (Pa)

dia = 1.15; %pipe inner diameter (m)
%rough = 0.05 * 1e-3; %steel (m)
rough = 0.0015 * 1e-3; %plastic (m)
lambda = 1 / (-2*log10(rough/(3.7 * dia)))^2; %darcy friction factor (-)
%lambda = 0.0;

area = dia.^2/4*pi; %(m2)
circum = dia * pi; %(m)
Kconcrete = 4; %thermal conductivity (W/(mK))
Ksteel = 50;

thick=4.3*2.54*0.01; %thickness of concrete wall
wattperkelvinpermeter = Kconcrete/thick*circum; % W/(mK)
wattperkelvinperm3 = wattperkelvinpermeter / area; %(W/(m^3K))

save_ = 1;
load_ = 0;

shocktube = 1;
model_temperature = 1; %model heat conduction from pipe walls (0 = adiabatic)
use_fric = 1; %include friction
choked_bc = 0; %force approximate choked flow at outlet

%duration of simulation (s)
endTime = 3600 * 24 * 3;
%endTime = 3600*5;
%endTime = 3600 * 10;
%endTime = 60 * 5;
%endTime = 60;
%endTime = 60/433;
%endTime = 400/433;

dtplot = endTime / 1000;
%dt_pplot = 1;
dt_pplot = endTime / 100;
CFL     = 0.3;
%% Grid Initialization

dx_init = 100; %length of control volume at opening (m)
dx_max = 1000; %maximum length of control volume (m)
fac = 1.1; %increase length of control volumes by this factor away from opening (-)

dx_c = -dx_init;
pos = 0;

x = zeros(1);

ns1_len = 1224 * 1000; %length of russian side of ns1(m)
ns2_len = 1230 * 1000;  %length of russian side of ns2(m)

ns1d = 230 * 1000; %length of german side of ns1 (m)
ns2d = 150 * 1000; %length of german side of ns2 (m)

pipe_length = (ns2_len - ns2d);
%pipe_length = (ns1_len - ns1d);
%pipe_length = (ns2_)

while (x(end) > -pipe_length)
    x = [x, x(end) + dx_c];
    if (dx_c > -dx_max)
        dx_c = dx_c * fac;
    end
end

x(end) = -pipe_length;



x = x(end:-1:1);

%x = linspace(0, 433 * endTime * 3,8000);

N = length(x) - 1;

%x = -x(end:-1:1) + x(end);
xc = 0.5 * (x(1:end-1) + x(2:end));

dx = diff(x);
time    = 0 ;
%% Initial Conidtions

if (shocktube)

    pressureR = p0_ambient;
    pressureL = p0_pipe;

    cfacR = methane_compression_factor(T0, pressureR);
    cfacL = methane_compression_factor(T0, pressureL);
    
    %densityR    = 0.1*0.72 ;       densityL    = 0.72 * 164 ;
    densityR    = pressureR / (T0 * Rgas) / cfacR;       densityL    = pressureL / (T0 * Rgas) / cfacL;
    
    cl = sqrt(gamma*pressureL/densityL);
    umid = 2 / (1 + gamma) * cl;
    
    fact = 1. - 0.5 * (gamma-1) * umid / cl;
        
    rho     = zeros(N,1) ;
    p       = zeros(N,1) ;
    u       = zeros(N,1) ;
    
    
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


T = ones(N,1) * T0;
e = T * Cv .* rho;

%%
new_rho = rho ;
new_u   = u   ;
new_e   = e   ;
new_T = T;
new_mom = rho.*u;

totmass = [];
mass_right = [];
mass_left = [];
e_right = [];
e_left = [];
plot_times = [];
tot_e = [];
p_base = [];
massflow_out = [];
u_out = [];
rho_out = [];
e_out = [];
T_out = [];
p_out = [];

pressure_plot = [];
e_plot = [];
u_plot = [];
rho_plot = [];
T_plot = [];
pplot_times = [];

plot_time = time;

pplot_time = time;

nstep = 0;

tstart = tic;

while time <= endTime
    

    nstep = nstep + 1;
    
    T = (e - 0.5*rho.*u.^2) ./ (Cv .* rho);
    %calculate compression factor based on previous time step pressure 
    %to avoid making this insanely complicated...
    
    for i = 1:1
        cfac = methane_compression_factor(T, p);
        p = T .* rho * Rgas .* cfac;
    end

    a = sqrt(cfac .* T *Rgas * gamma); %speed of sound
    
    dt      = min(CFL./(a(:)+abs(u(:))) .* dx(:)) ;  % cfl condition time step
    time    = time + dt ;

    if (1)
        mom = rho .* u;
        lambda_t = abs(u) + a;

        
        lambda_e = max(lambda_t(1:end-1), lambda_t(2:end)); %numerical diffusion

        mom_flux_t = rho.*u.*u + p;
        energy_flux_t = u .* (e + p);

        rho_flux_e = 0.5 * (mom(1:end-1) + mom(2:end)) + 0.5 * lambda_e .* (rho(1:end-1) - rho(2:end));
        mom_flux_e = 0.5 * (mom_flux_t(1:end-1) + mom_flux_t(2:end)) + 0.5 * lambda_e .* (mom(1:end-1) - mom(2:end));
        energy_flux_e = 0.5 * (energy_flux_t(1:end-1) + energy_flux_t(2:end)) + 0.5 * lambda_e .* (e(1:end-1) - e(2:end));

        

        if (choked_bc)
            if (u(N-1) < a(N-1) && u(N - 1) > 0)
                vc = u(N - 1) / a(N - 1);

                xd = vc * (1 + gamma) / (2 + vc * (gamma - 1));

                uoutlet = u(N - 1) / xd;

                c0 = 0.5 * (1 + gamma) * uoutlet;

                p0 = p(N - 1) * (c0 / a(N - 1)) ^ (2*gamma/ (gamma - 1));

                poutlet = p0 * (2 / (1 + gamma))^(2*gamma/(gamma - 1));
                rhooutlet = gamma * poutlet / uoutlet^2;
                eoutlet = uoutlet.^2*rhooutlet/2 + poutlet / (gamma - 1) / cfac(N - 1);

                Toutlet = (eoutlet - 0.5*rhooutlet.*uoutlet.^2) ./ (Cv .* rhooutlet);

                if (poutlet > pressureR)
                    rho_flux_e(end) = rhooutlet * uoutlet;
                    energy_flux_e(end) = uoutlet * (eoutlet + poutlet);
                    mom_flux_e(end) = rhooutlet * uoutlet^2 + poutlet;
                end

            elseif (u(N - 1) == 0)
                uoutlet = 2 / (1 + gamma) * a(N - 1);
                poutlet = p(N - 1)*(2 / (1 + gamma)) ^ (2 * gamma / (gamma - 1));
                rhooutlet = gamma * poutlet / uoutlet^2;

                eoutlet = uoutlet.^2*rhooutlet/2 + poutlet / (gamma - 1) ./ cfac(N-1);

                Toutlet = (eoutlet - 0.5*rhooutlet.*uoutlet.^2) ./ (Cv .* rhooutlet);

                rho_flux_e(end) = rhooutlet * uoutlet;
                energy_flux_e(end) = uoutlet * (eoutlet + poutlet);
                mom_flux_e(end) = rhooutlet * uoutlet^2 + poutlet;
            end
        end


        new_rho(2:end-1) = rho(2:end-1) + dt ./ dx(2:end-1)' .* (rho_flux_e(1:end-1) - rho_flux_e(2:end));
        new_mom(2:end-1) = mom(2:end-1) + dt ./ dx(2:end-1)' .* (mom_flux_e(1:end-1) - mom_flux_e(2:end));
        new_e(2:end-1) = e(2:end-1) + dt./dx(2:end-1)' .* (energy_flux_e(1:end-1) - energy_flux_e(2:end));


        %closed BC on left side of pipe
        new_rho(1) = rho(1) - dt/dx(1) * rho_flux_e(1);
        mom1 = rho(1)*u(1) - dt/dx(1) * (mom_flux_e(1) - (0*rho(1)*u(1).^2 +p(1)));
        new_u(1) = mom1 / new_rho(1);
        new_e(1) = e(1) - dt/dx(1) * (energy_flux_e(1) - 0*(u(1) * ( e(1) + p(1) )));

        new_u(2:end-1) = new_mom(2:end-1) ./ new_rho(2:end-1);

        %u_fric = (2 * new_u) ./ (1+sqrt(1+(16*dt.*abs(new_u) * lambda) / dia));
        
        if (use_fric)
            %implicit formula for friction to ensure stability
            
            u_fric = (2 * new_u) ./ (1+sqrt(1+(2*dt.*abs(new_u) * lambda) / dia));
            new_u = u_fric;
            
            %explicit friction
            %fric = new_rho .* sign(new_u).*new_u.^2 / 8 * lambda ./ dia * 4.; %dimension Pa/m
            %fric_mom = new_u .* new_rho - dt * fric;
            %new_u = fric_mom ./ new_rho;
        end

        
        


    end
    

    if (model_temperature)
        new_T = (new_e - 0.5.*u.*u.*rho) ./ (rho * Cv);
        tdiff = T0 - new_T;

        deltaT = tdiff .* (1 - exp(-dt * wattperkelvinperm3 ./ (rho*Cv)));
        deltaE = deltaT .* rho .* Cv;

        new_e = new_e + deltaE;
    end

    %approximate open BC
    new_u(N) = new_u(N - 1);
    new_rho(N) = new_rho(N - 1);
    %force pressure to ambient pressure in last cell
    new_e(N) = pressureR / (gamma-1) / cfac(end) + 0.5*new_rho(N) .* new_u(N) .* new_u(N);


    %save variables
    if (time >= plot_time)
        fprintf("time: %f\n", time);
        telapsed = toc(tstart);
        fprintf("steps per second: %f\n", nstep / telapsed);
        plot_time = plot_time + dtplot;
        plot_times = [plot_times time];
        totmass = [totmass sum(rho(:).*dx(:).*area)];
        mass_right = [mass_right sum(rho(floor(N/2) + 1:end).*dx(floor(N/2) + 1:end)'.*area)];
        mass_left = [mass_left sum(rho(1:floor(N/2)-1).*dx(1:floor(N/2)-1)'.*area)];
        tot_e = [tot_e sum(e)];
        e_left = [e_left sum(e(1:floor(N/2)-1))];
        e_right = [e_right sum(e(floor(N/2):end))];
        p_base = [p_base p(1)];
        %massflow_out = [massflow_out rho_flux_e(end)];
        massflow_out = [massflow_out u(N - 1) * rho(N - 1)];
        u_out = [u_out u(N - 1)];
        rho_out = [rho_out, rho(N - 1)];
        e_out = [e_out e(N - 1)];
        T_out = [T_out T(N - 1)];
        p_out = [p_out p(N - 1)];
    end

    if (time >= pplot_time) 
        pplot_time = pplot_time + dt_pplot;
        pplot_times = [pplot_times time];
        pressure_plot = [pressure_plot p];
        e_plot = [e_plot e];
        u_plot = [u_plot u];
        rho_plot = [rho_plot rho];
        T_plot = [T_plot T];
    end
    
     rho     = new_rho ;
     u       = new_u ;
     e       = new_e ;
     
     
end

if (save_)
    save("workspace.mat");
    rho_plot_save = rho_plot;
    e_plot_save = e_plot;
    T_plot_save = T_plot;
    pressure_plot_save = pressure_plot;
    u_plot_save = u_plot;
end

