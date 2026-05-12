%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           Reaction Rate_C                              %%
%% written by: ems                                                       %%
%% Feb 13, 2025                                                          %%
%% -calculates the rate of C production at a given P-T condition         %%
%% -note: this assumes you already have already found (i,j) of the P-T   %%
%%  conditions. Called by degas2.m                                        %%
%%                                                                       %%
%% -rxn is:                       gph + H2O =                            %%
%%                            0.5 CH4 + 0.4 CO2                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dmC=Rxnrate_C(T,GtotC,i,j)

%output dmC is in kg of C per m3 rock per year

%input variables: T in Celsius
                % Gtot is DGrxn gridded to P-T range. Calculated offline
                        % in theriak-domino.
                % (i,j) are the indices to call Gtot at the right P-T 

%define constants
A=1;            % surface area per volume, cm2 cm-3; after Ague & Rye 1998
TK=T+273.15;    %convert T input to Kelvin
nu=1;          %stoichiometric coefficient of (CO2+CH4) in reaction 
R=8.314;        %gas constant J mol^-1 K-1

    %________turn on for Ague and Rye 1998_______%
    order=2.68;     % reaction order
    Ea=83700;       % Activation Energy after Ague & Rye 1998
    ko=1.38*10^-14; % rate constant after Ague & Rye 1998
    Trefinv=1/873.15; % Reference Temperature in Kelvin, after Ague & Rye 1998

    %________turn on for Lasaga and Rye 1993_______%
    %order=1;        % reaction order
    %Ea=83700;       % Activation Energy after Lasaga & Rye 1993  
    %ko=1.051*10^-8; % rate constant after Lasaga & Rye 1993 
    %Trefinv=1/873.15; % Reference Temperature in Kelvin, after Lasaga & Rye 1993
    
    %________turn on for Chang and Berner 1999_______%
    %order=0.5;        % reaction order
    %Ea=43000;         % Activation Energy after Chang & Berner 1999
    %ko=6.3072*10^-9;  % rate constant after Chang & Berner 1999
    %Trefinv=1/(24+273.15);   % Reference Temperature after Chang & Berner 1999
    
    %________turn on for Burnham & Sweeney 1990_______%
    %order=0;           % reaction order
    %Ea=221752;          % Activation Energy after Burnham & Sweeney 1990
    %ko=1.313*10^18;    % rate constant after Burnham Sweeney
    %Trefinv=0;         % Reference temperature of infinity after Burnham Sweeney
    %nu=0.425;          % update nu for 50% conversion of 0.85 
    

%% find DGrxn
DG=GtotC(i,j);   % find DG of reaction at (i,j)

%% calculate reaction rate
k=ko*exp(-Ea/R*(1/TK-Trefinv));  % Arrhenius equation after Ague & Rye 1998

dmdt=k*nu*A*(abs(DG))^order;    % rxn rate in moles C per cm3 per year

dmdt=dmdt*100*100*100;          % convert to moles C m-3 yr-1

dmdt=dmdt*12.01;                % convert to g C m-3 yr-1

dmC=dmdt/1000;                  % convert to kg C m-3 rock yr-1

end
