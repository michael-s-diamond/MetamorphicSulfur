%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Thermal Modeling of Sill Intrusion                  %%
%% written by: ems                                                      %%
%% date: Jan 24, 2022                                                   %%
%% updated for Sulfur Feb 21, 2025                                      %%
%% calculates the changing temperature as a result of sill intrusion    %%
%% using a finite difference second derivative scheme                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,dtyrs,z,dz,x,Tdiff]=thermaldiffusion()

% output:     time is vector in years, dtyrs timestep in years
            % z is depth vector in m, dz depth step in m
            % T is matrix of temperatures at each time and z in Celsius
            % x is the rock type (0= country rock, other value = sill)
%% set up depth variables 
z1=0;               %minimum depth (m)
z2=4000;            %maximum depth (m)
dz=1.5;               %depth resolution (m)
z=[z1:dz:z2]';      %matrix to hold all depths
zsize=size(z,1);
    
%% set rock types; 0=country rock; anything else is a sill (numbered 1 +)
x=zeros(zsize,1);                       % holds rock type per depth

for i=1:zsize                           % find sills and number 1 to n
    if  z(i)>260 && z(i)<(260+143)      % this is manually hard-coded to
        x(i)=1;                         % to match the stratigraphy of
    elseif z(i)>533 && z(i)<(533+48)    % borehole OAST after Heimdal 2018
        x(i)=2;
    elseif z(i)>1408 && z(i)<1408+4 
        x(i)=3;
    elseif z(i)>2177 && z(i)<(2177+3)
        x(i)=4;
    elseif z(i)>2784 && z(i)<(2784+62) 
        x(i)=5;
    end
end

%% set intrusion timing 
intrusion_ages=[201.5*10^6, 201.499*10^6, 201.498*10^6,201.497*10^6, ...
    201.496*10^6];      % sill intrusion ages in Ma,every 1000 years

flowtime=10;                    %sill duration in years
flowtime=flowtime*525600*60;    % sill duration in seconds

timecounter=flowtime+1;

%% set initial Temperature
%  initial geotherm from OAST (Heimdal 2018) is 70C at 2000m and 0 at surf.
                        
dTdz=70/2000;               % geothermal gradient in  degrees C per meter
To=z*dTdz;                  % starting temperature per depth via geotherm

Tigno=1100;                 % starting temperature of igneous intrusion


%% set diffusion constant in m^2 s^-1
K=10^-6;

%% put time values in seconds for diffusion calculation 

dt=1000000;                         % set time step in seconds
dtyrs=dt/(3.154*10^7);              % convert timestep to yrs for later
    
tstartyears=201.501*10^6;           % start time in years ago
tstartsec=tstartyears*3.154*10^7;   % start time in seconds ago 
tendyears=201.493*10^6;             % end time in years ago
tendsec=tendyears*3.154*10^7;       % end time in seconds ago

duration=tstartsec-tendsec;         % total duration in seconds
tsize=ceil(duration/dt)+1;          % calculate total number of timesteps

time=[tstartyears:-dtyrs:tendyears]; %hold the age of each timestep in yrs

Tdiff=zeros(zsize,tsize);        % hold temperature profile info at each time
Tdiff(:,1)=To;                   % use initial T profile at start time

%% make beta term and test for stability
B=K*dt/(dz^2);

if B>=0.5 
    fprintf('unstable: check step sizes');
end


%% Main Thermal Diffusion Loop

for j=2:1:tsize             % loop over time
    
    prev=Tdiff(:,j-1);          % hold T profile at previous timestep
    
    %check for new sill intrusion
    
    
    [k,l]=find(intrusion_ages==time(j)); % check if age = age of intrusion
     if ~isempty(k)                      % if yes
         [n,m]=find(x==l);               % find  depth range of sill # l
         prev(n,1)=Tigno;                % set T=T igneous within sill
         timecounter=0;
         Lperm=l;
     end   
   
     if timecounter <= flowtime
         [n,m]=find(x==Lperm);
         prev(n,1)=Tigno;
     end
         
     
    for i=2:(zsize-1)       % loop over depth        
        Tdiff(i,j)=prev(i+1)*B+prev(i)*(1-2*B)+prev(i-1)*B; % 2nd deriv. scheme
    end
    
    %set boundary condition:top is Dirichlet, fixed at initial value
    Tdiff(1,j)=prev(1);
    
    %bottom is Dirichlet, fixed at initial value 
    Tdiff(zsize,j)=prev(zsize);
    
    timecounter=timecounter+dt;
end

%% change x so it only holds a value if the rock can degas S
for i=1:zsize
    if x(i)==0 
            x(i)=1;
    else
            x(i)=0;
    end
end

%%  make figures 
figure
plot(z,Tdiff(:,j)); % plot final T profile
hold
xlabel('Temperature Celsius');
ylabel('depth (m)')
Tmax=max(Tdiff');   
plot(z,Tmax);   % plot maximum T profile
end


