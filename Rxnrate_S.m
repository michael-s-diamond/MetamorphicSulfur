%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           Reaction Rate_S                             %%
%% written by: ems                                                       %%
%% Feb 13, 2025                                                          %%
%% -calculates the rate of S production at a given P-T condition         %%
%% -note: this assumes you already have already found (i,j) of the P-T   %%
%%  conditions. Called by degas2.m                                       %%
%%                                                                       %%
%% -rxn is:              Pyr + 0.5 gph + H2O =>                          %%
%%          0.91 H2S + 0.5 CO2 +0.1 H2 + 0.4408 trot + 0.7008 trov       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dmS=Rxnrate_S(T,Gtot,order,i,j)

%output dmS is in kg of S per m3 rock per year

%input variables: T in Celsius
                % Gtot is DGrxn gridded to P-T range. Calculated offline
                        % in theriak-domino.
                % order is the reaction order and must be either 1 or 2.68
                % (i,j) are the indices to call Gtot at the right P-T 

%define constants
TK=T+273.15;    %convert T input to Kelvin
nuS=0.91;       %stoichiometric coefficient of H2S in reaction 
R=8.314;        %gas constant J mol^-1 K-1
                   
A=1;            % surface area per volume, cm2 cm-3; after Ague & Rye 1998
 
%% find DGrxn
DG=Gtot(i,j);   % find DG of reaction at (i,j)

if DG > 0       % if the reaction would go backwards
    dmS=0;      % set rate = 0
    return      % and return control to degas.m     
end

%% establish reaction order
if order == 1
    Ea=83700;       % Activation Energy after Lasaga & Rye 1993  
    ko=1.051*10^-8; % rate constant after Lasaga & Rye 1993 
    Trefinv=1/873.15; % Reference Temperature in Kelvin, after Lasaga & Rye 1993

                            % units of mols cm-2 yr-1 (J/mol)-1
elseif order == 2.68
    Ea=83700;       %Activation Energy after Ague & Rye 1998
    ko=1.38*10^-14;  % rate constant after Ague & Rye 1998
    Trefinv=1/873.15; % Reference Temperature in Kelvin, after Ague & Rye 1998

                            % units of mols cm-2 yr (J/mol)-2.68
elseif order == 0
    %________turn on for Hong and Fegley 1997_______%
    Ea=297000;       % Activation Energy after Hong and Fegley 1997 (J/mol)
    ko=1.01*10^20;  % rate constant after Hong and Fegley 1997, mols cm-2 yr
    Trefinv=0;      % Reference temperature of infinity after Hong and Fegley 1997 
    
    %_____turn on for Fegley et al._______%
    %Ea=156000;      % Activation Energy after Fegley et al (J/mol)
    %ko=0.011679;    % rate constant after Fegley et al (mols cm-2 yr)
    %Trefinv=1/663;   % Reference temperature after Fegley et al (K)
    
     %_____turn on for Lambert et al. in inert gas_______%
    %Ea=297000;      % Activation Energy after Lambert et al (J/mol)
    %ko=1.168*10^18; ;    % rate constant after Lambert et al (mols cm-2 yr)
    %Trefinv=0;         % Reference temperature of infnity after Lambert et al (K)
    
     %_____turn on for Lambert et al. in H2 at 0.1 Mpa _______%
    %if TK < 670
    %    ko=4.43475;     % rate constant after Lambert et al (mols cm-2 yr)
    %    Ea=31000;       % Activation Energy after Lambert et al (J/mol)
   %elseif TK < 780
   %   ko=139612;       % rate constant after Lambert et al (mols cm-2 yr)
   %   Ea=89000;        % Activation Energy after Lambert et al (J/mol)
  %else 
  %    ko=1.807*10^12;  % rate constant after Lambert et al (mols cm-2 yr)
  %    Ea=195000;       % Activation Energy after Lambert et al (J/mol)
  %end  
  %Trefinv=0;           % Reference temperature of infinity after Lambert et al (K)
   
    
else
    fprintf('uh oh! this reaction order not supported! pick 0, 1 or 2.68\n');
    return
end

%% calculate reaction rate
k=ko*exp(-Ea/R*(1/TK-Trefinv));  % Arrhenius equation after Ague & Rye 1998

dmdt=k*nuS*A*(abs(DG))^order;   % rxn rate in moles H2S per cm3 per year

dmdt=dmdt*100*100*100;          % convert to moles H2S (or S) m-3 yr-1

dmdt=dmdt*34;                   % convert to g S m-3 yr-1

dmS=dmdt/1000;                  % convert to kg S m-3 rock yr-1

end
