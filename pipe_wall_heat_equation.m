Ny = 100;     %number of radial temp measurement points
Nt = 110;     %number of time steps
density = 7700;  %kg/m3 of stainless steel
heat_capacity = 450;   %J/kg/K
thermalconductivity =14;  %J/K/m/s

TimePeriodSeconds = 600;  %total time including warm time
SteelThickness = 0.027;      %meters
ConcreteThickness = 4.3*2.54*0.01; %meters
Depth = SteelThickness + ConcreteThickness;    %thickness of pipeline steel (meters)

Tenv = 285;
Tgas = 230;

cSteel = thermalconductivity/(heat_capacity*density) %diffusivity m^2/s
cConcrete =  9.51e-7 %thermal diffusivity of concrete

dt = TimePeriodSeconds / (Nt - 1);
dx = Depth / (Ny + 1);

WaterTemp = 12;
GasTemp   = -50;
simulationTime=600;  %max duration of cold  

%[~,~,A] = laplacian(Ny,{'DD'});

diffusion=ones(1,Ny + 1);
iSCboundary= floor(Ny*SteelThickness/Depth);

diffusion(1:iSCboundary - 1)=cSteel;
diffusion(iSCboundary) = 2. / (1 / cSteel + 1 / cConcrete); %harmonic mean
diffusion(iSCboundary + 1:end) = cConcrete;


diffusion=(dt/dx^2)*diffusion;

Ads = spdiags([-diffusion(2:end)' (diffusion(1:end-1)' + diffusion(2:end)') -diffusion(1:end - 1)'], [-1, 0, 1], Ny, Ny);
Ad = (Ads + speye(Ny));


e0=zeros(Ny, 1);
e0(1)=Tgas;

e1 = zeros(Ny, 1);
e1(end) = Tenv;

T = zeros(Ny, Nt);
T(:, 1) = ones(Ny, 1) * Tenv;

for i=1:Nt-1
         T(:, i+1) = Ad\(T(:,i) - Ad(2,1) * e0 - Ad(end-1,end) * e1);
end

millimeters = [0. linspace(0.5*dx, Depth - 0.5*dx, Ny) Depth];

for seconds=0:60:simulationTime
    
    step  = floor(seconds / dt) + 1;

    Tplot = [Tgas; T(:, step); Tenv] - 273;

    plot(millimeters, Tplot);

    title('Steel temperature every minute first 10 minutes')
    xlabel('millimeter from inner surface')
    ylabel('Degree Celsius temp')
    hold on
end
% for seconds=0:60:simulationTime
% 
%     BoundaryTemp = zeros(Nt,1);
%     pulses = [0,seconds,(GasTemp -WaterTemp)]';       %Cooling pulse
% 
%     pulses= [pulses, [seconds+1,simulationTime,0]'];         % Before explosion +12C
% 
% 
% 
%     for i=1:size(pulses,2)
%         pu=pulses(:,i);
%         iPs = floor(Nt* pu(1)/TimePeriodSeconds)+1;
%         iPe = floor(Nt* pu(2)/TimePeriodSeconds);
%         BoundaryTemp(iPs:iPe) = pu(3);
%     end
% 
% 
%     Dtemp = B*BoundaryTemp;  % cross section temp profile resulting from the pulses
% 
%     increment = 1000*Depth/(Ny-1);
%     millimeters = 0: increment : (1000*Depth); %horizontal scale in plot
% 
%     plot(millimeters, Dtemp+WaterTemp);
% 
% 
%     title('Steel temperature every minute first 10 minutes')
%     xlabel('millimeter from inner surface')
%     ylabel('Degree Celsius temp')
%     hold on
% end
%   grid on
% 
% hold off

