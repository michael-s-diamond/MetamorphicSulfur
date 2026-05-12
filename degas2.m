%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   Degas2                              %%
%% written by: ems                                                       %%
%% Feb 24, 2025                                                          %%
%% calculate the Sulfur degassed given input of T at each depth and time %%
%% note C degassed is just forced to match S given rxn stoichiometry     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sflux,Cflux,massSz,massCz]=degas2(Pm,Tm,massSrem,massCrem,...
                                            z,Tdiff,time,Gtot,GtotC,x)

%outputs      Sflux is outgassing of S in kg per year
            % Cflux is outgassing of CO2 (not C) in kg per year
            % massSz: mass of S remaining at each depth and time in kg
            % massCz: mass of CO2 remaining at each depth and time in kg
            
%inputs       Pm and Tm: P and T meshgrid used to call thermo data
            % Gtot is the DG of S releasing rxn at each P-T in J
            % GtotC is the DG of C releasing rxn at each P-T in J
            % massSrem: mass of S remaining at equilibrium at each P-T
                                                     %units of kg m-3
            % massCrem: mass of C rem at eq at each P-T; units of kg m-3                                         
            % z: depth data from diffusion simulation in m
            % Tdiff: Temp. in C at each z and time from diffusion simulation
            % time: time data from diffusion simulation in years ago
            % dt: size of timestep in diffusion simulation
            % x: rock types in diffusion simulation
            

%% first lower the resolution on Tdiff to make calculation manageable

dtnew=1;    % new size of timestep in years
dznew=1;    % new size of depth step in m

znew=[min(z):dznew:max(z)]';            %z data at new resolution
tnew=[max(time):-dtnew:min(time)];      %t data at new resolution

Tdiff=interp2(time,z,Tdiff,tnew,znew);  %interpolate T to new resolution
x=interp1(z,x,znew,'nearest');          %interpolate x to new resolution 
 
%% initialize 

zsize=size(Tdiff,1);                    %holds the size of z vector
tsize=size(Tdiff,2);                    %holds the size of t vector

massSz=zeros(size(Tdiff));              %holds mass of S (see def. above)     
massCz=zeros(size(Tdiff));              %holds mass of C (see def. above)


massSz(:,1)=max(max(massSrem))*x; %set S content at t=0
massCz(:,1)=max(max(massCrem))*x; %set C content at t=0

order=0;                          %define S reaction order. for S Must be 0, 1 or 2.68


%% main time loop
for j=2:1:tsize                   %loop over time
    for i=1:1:zsize               %loop over depth
        
        P=1+round(znew(i)*0.098); %Calculate hydrostatic 98 bars/km
        T=round(Tdiff(i,j));      %find T at timestep j and depth i
        
%% start S calculation
        if T<200                  %if T < 200 no need to do calculatio
            massSz(i,j)=massSz(i,j-1); %no mass change from previous t 
            massCz(i,j)=massCz(i,j-1); %no change from previous t
        else
            %[n,m]=find(Pm==P & Tm==T);% find indices of P and T values
            %replacing the above with below for efficiency
            n=P;
            m=T-199;
            
            Seq=massSrem(n,m);        % find S remaining at this P-T at eqm
      
            dmS=Rxnrate_S(T,Gtot,order,n,m); %calculate S rxn rate at this T
            massSz(i,j)=massSz(i,j-1)-dmS*dtnew;  %find new mass of S
        
            if massSz(i,j) < Seq          % check for passing the Seq #
               massSz(i,j)=Seq;           % and reset to Seq # if passed
            end
        
            if massSz(i,j)>massSz(i,j-1)   % check for re-growing S in rock
               massSz(i,j)=massSz(i,j-1);  % and reset to previous value
            end


%% start C calculation
        Ceq=massCrem(n,m);          %find C remaining at this P-T at eqm
        dmC=Rxnrate_C(T,GtotC,n,m); %calculate C rxn rate at this T
        
        massCz(i,j)=massCz(i,j-1)-dmC*dtnew;  %find new mass of C
        
         if massCz(i,j) < Ceq          % check for passing the Ceq #
               massCz(i,j)=Ceq;           % and reset to Ceq # if passed
            end
        
            if massCz(i,j)>massCz(i,j-1)   % check for re-growing C in rock
               massCz(i,j)=massCz(i,j-1);  % and reset to previous value
            end
        
        end                         % end if T < 200        
    end                             % end loop over depth 
 
end                                 % end loop over time

%% calculate flux to surface

masssumS=sum(massSz)*dznew;    %integrate rock column for mass of S (kg/m2)
masssumC=sum(massCz)*dznew;    %integrate rock column for mass of C (kg/m2)

A=10^5;                        % total area in km2
A=A*1000*1000;                 % total  are in m2

Sflux=A*(-diff(masssumS)/dtnew);    % mass flux to the surface in kg S/yr   
Cflux=A*(-diff(masssumC)/dtnew);    % mass flux to the surface in kg C/yr

end                                 % end function
    
