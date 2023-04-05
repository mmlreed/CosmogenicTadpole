% Script for computing the coevolution of topography and soil chemical
% weathering, associated with the manuscript in review, "The importance of
% hillslope scale in responses of chemical erosion rate to change in
% tectonics and climate", by Ken Ferrier and Taylor Perron.
% 
% Additional cosmogenic nuclide functionality, nonlinear hillslope transport, and ability to work 
% with more complex channel networks by Miles Reed for the paper "

% Input grids and parameters:
%
%   B: Matrix of initial bedrock surface elevation
%   H: Matrix of initial soil thickness
%   X: K x J x nX array of initial mineral concentrations (where n is the
%      number of mineral species)
%   C: Matrix specifying location of 'boundary' cells where soil erosion, 
%      soil production, and weathering are zero
%   run_name: Character string naming the run. If specified, the parameters 
%             and elevations at each timestep will be saved in a binary 
%             file called <run_name>.mat
%   p: A struct array of parameters that includes:
% 
%     p.dt             Timestep (yr)
%     p.tf             Total time of the simulation (yr)
%     p.dx, p.dy       Grid spacing in the x and y directions (m)
%     p.D              Hillslope diffusivity (m^2/yr)
%     p.K              Coefficient in stream power-type incision law 
%                       (kg m^(1+2m) yr^-2)
%     p.tauc           Threshold for fluvial incision
%     p.m, p.w         Exponents on drainage area and slope, respectively
%     p.Kw             Coefficient in channel width power law: 
%                       width = Kw*A^wexp
%     p.wexp           Exponent in channel width power law
%     p.E              Rate of relative uplift or base level lowering (m/yr)
%     p.saveint        Grids are saved every saveint iterations. 
%                      No output will be saved if save_interval==0.
%     p.plotint        Plot will be redrawn every plotint iterations
%     p.zrange         Range of vertical axis of plots
%     p.vexag          Vertical exaggeration of output plots
%     p.P              Soil production rate at zero soil depth 
%                       (kg m^-2 yr^-1)
%     p.alpha          Decay constant of soil production rate with soil 
%                       depth (m^-1)
%     p.rhor           Rock bulk density (kg m^-3)
%     p.rhos           Soil bulk density (kg m^-3)
%     p.wx             Mineral molar mass (kg mol^-1)
%     p.kx             Mineral reactivity (mol m^-2 yr^-1)
%     p.sx             Clay production rate (mol m^-3 yr^-1)
%     p.Ax             Mineral specific surface area (m^2 mol^-1)
%     p.xr             nX-length vector of mineral concentrations in 
%                       bedrock (dimensionless)
%     p.xname          Cell array of mineral names, e.g. {'quartz,
%                       'plagioclase','K-feldspar'}
%     p.x2plot         Index of mineral to plot (1 <= x2plot <= nX)
%     p.Hchannelvalue  Soil thickness assigned to channels (m).  A very
%                       small value greater than zero is chosen here to 
%                       avoid divide-by-zero issues.
%     p.Estepfactor    Factor by which rock uplift rate increases in
%                       experiments with imposed changes in rock uplift
%                       rate
%     p.Pstepfactor    Factor by which soil production rate increases in
%                       experiments with imposed changes in soil production
%                       rate
%     p.Dstepfactor    Factor by which soil transport rate increases in
%                       experiments with imposed changes in soil production
%                       rate
%     p.kxstepfactor   Factor by which dissolution rates increase in
%                       experiments with imposed changes in dissolution
%                       rate
%     p.mark           criterion for approaching new steady state after
%                       step
%     p.n_atstep       Time step at which to impose perturbation
%     p.experiment_type   String indicating the type of imposed
%                          perturbation.  Must be one of 'tectonic',
%                          'climatic', 'D', 'P', 'kx', or 'steady_state'.
%


%% Load input grids and parameter values.
% Manually load grids from .mat -- Thi

%% Create perturbation.

% Here is an example for a step perturbation in rock uplift rate.
p.Estepfactor = 2;  % factor by which this changes in the imposed step.
% E.g., a factor of 1.5 would be a 50% increase.
p.experiment_type = 'tectonic';
p.tf = 5e6;  % duration of simulation (years)
run_name = 'tectonic';

% Save inputs as struct for the output.
% inputs.B = B;
% inputs.C = C;
% inputs.H = H;
% inputs.X = X;
% % Once we have an intial run from steady state, uncomment these.
% inputs.Nzb = Nzb;
% inputs.Ns = Ns;
% inputs.N_prof = N_prof; % steady state conc. profile needs inputed as it will be
% inputs.N_prof2 = N_prof2;
% inputs.p = p;
% inputs.run_name = run_name;

% Save outputs for do_output = 1.  Don't save outputs for do_output = 0.
do_output = 1;

%% PARAMETERS
% Parse the parameters input vector
dx = p.dx; 
dy = p.dy;
D = p.D;
E = p.E;
P = p.P;
alpha = p.alpha;
rhor = p.rhor;
rhos = p.rhos;
kx = p.kx;
Ax = p.Ax;
sx = p.sx;
wx = p.wx;
xr = p.xr;
host_min = p.host_min;
zr_min = p.zr_min;
lambda_10Be = p.lambda_10Be;
L1 = p.L1;
L2 = p.L2;
L4 = p.L4;
xname = p.xname;
x2plot = p.x2plot;
mark = p.mark;  % criterion for approaching new steady state after step
experiment_type = p.experiment_type;

% Define x-y grids and save as parts of inputs.
[x, y] = meshgrid(0:dx:dx*(size(C,1)-1), 0:dy:dy*(size(C,1)-1));
inputs.x = x;
inputs.y = y;

% Don't save output if p.saveint <= 0.
if round(p.saveint) <= 0
    do_output = 0;
else
    save_interval = round(p.saveint);
end

%% TIME VARIABLES
t = 0; % time in yr; this is the time at the first timestep
dt = 10; % delta t
N = round(p.tf/dt); % number of iterations to do                           
dt_initial = p.dt;  % save this in case you change dt mid-simulation

%% INITIAL CONDITIONS
[K, J, nX] = size(X); % nX is the number of mineral species

%% Conditions for imposed step change.

switch experiment_type
    
    case 'tectonic'  % consider a perturbation in channel lowering rate
        E_afterstep = E * p.Estepfactor;
        
    case 'climatic'  % consider simultaneous perturbations in D, P, kX
        K3_afterstep = K3 * 1.5;
        P_afterstep = P * 1.5;
        kx_afterstep = kx * 1.5; % 150% increase in these params
        
    case 'steady_state'  % change nothing
        
    otherwise
        error('Choose one of the valid experiment types')

end

%% MEMORY ALLOCATION

Xnminus1 = zeros(K,J,nX); % value of X at time n-1 (used in leapfrog scheme)
Hnminus1 = zeros(K,J); % used in Adams-Bashforth
Bnminus1 = zeros(K,J); % used in Adams-Bashforth

Xtemp = zeros(K,J,nX); % these three are used to keep track of X,H,B at time n-1
Htemp = zeros(K,J); 
Btemp = zeros(K,J); 

% Assign values to first timestep
if do_output % if run_name was specified, initialize grids
    
	output.topo_m = zeros(K,J,ceil(N/save_interval)+1); % each 'layer' of
        % this data cube holds the matrix of elevations at time n.
	output.H_m = zeros(K,J,ceil(N/save_interval)+1);  % soil thickness
	output.X = zeros(K,J,nX,ceil(N/save_interval)+1);  % mineral abundances
    output.W_mperyr = zeros(K,J,ceil(N/save_interval)+1);  % chemical erosion rate
    output.R_mperyr = zeros(K,J,ceil(N/save_interval)+1);  % physical erosion rate
    output.Dact_mperyr = zeros(K,J,ceil(N/save_interval)+1); % both phys and chem 
    output.Dinf_mperyr = zeros(K,J,ceil(N/save_interval)+1);
    output.Dinf_lal_mperyr = zeros(K,J,ceil(N/save_interval)+1); % inferred d rate from cosmo
    output.Dinf_zb_mperyr = zeros(K,J,ceil(N/save_interval)+1);
    output.Dalt_mperyr = zeros(K,J,ceil(N/save_interval)+1);
    output.Winf_mperyr = zeros(K,J,ceil(N/save_interval)+1);
    output.Einf_mperyr = zeros(K,J,ceil(N/save_interval)+1);
    output.Winf_alt_mperyr = zeros(K,J,ceil(N/save_interval)+1);
    output.Einf_alt_mperyr = zeros(K,J,ceil(N/save_interval)+1);
    output.Ns_atg = zeros(K,J,ceil(N/save_interval)+1);
    output.Nzb_atg = zeros(K,J,ceil(N/save_interval)+1);
    
    % Save values at first save timestep.
    output.topo_m(:,:,1) = B + H;  % topography (m)
    output.H_m(:,:,1) = H;  % soil thickness (m)
    output.X(:,:,:,1) = X;  % soil mineral abundances ()
    output.t_yr(1) = t;  % time since onset of simulation (yr)
    output.Ns_atg(:,:,1) = Ns/1000; % Ns in atoms g-1;
    output.Nzb_atg(:,:,1) = Nzb/1000; 
    output.Dact_mperyr(:,:,1) = Dact;
    output.R_mperyr(:,:,1) = R;
    output.W_mperyr(:,:,1) = W;
    output.Dinf_mperyr(:,:,1) = Dinf;
    output.Dinf_lal_mperyr(:,:,1) = Dinf_lal;
    output.Dinf_zb_mperyr(:,:,1) = Dinf_zb;
    output.Einf_mperyr(:,:,1) = Einf;
    output.Winf_mperyr(:,:,1) = Winf;
    output.Dalt_mperyr(:,:,1) = Dalt; % alternative 
end


%% MAIN ITERATION LOOP

n_atstep = p.n_atstep;  % timestep at which step occurs
tstep = dt * n_atstep;  % time (yr) at the step iteration
inputs.timeatstep_yr = tstep;  % save time of step (yr)

for n=1:N  % timesteps (N = total number of timesteps)
    
    % This iterative loop computes changes in H, B, and X in a splitting
    % method, which computes portions of the changes in each term in
    % sequential steps rather than all at once. For example, part of the 
    % change in H from one timestep to the next is computed in one function
    % (SoilProd), and the updated value of H is fed into a second function
    % (Weather), which updates H again and feeds the newly updated H into a
    % third function (Erode).  Such splitting methods are favored here
    % because of their higher accuracy over non-splitting methods.
    
    % Impose step perturbation.  Note that if you change D here, this is
    % also where you'll need to reassign values for Aplus, Aminus, Bplus,
    % and Bminus too.
    
    if n == n_atstep
        
        switch experiment_type
            
            case 'tectonic'  % change only channel lowering rate E
                E = E_afterstep;
            
            case 'climatic'  % simultaneously change D, P, and kX
                P = P_afterstep;
                kx = kx_afterstep;
                K3 = K3_afterstep;  % requires changing Aplus Aminus Bplus
            case 'steady_state'  % change nothing
        end
        
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = t + dt;
    Htemp = H;
    Btemp = B;
    Xtemp = X;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%% Soil Production %%%%%%%%%%%%%%%%%%%%%%%%
    if n == 1
        BLV_minus1 = BLV; 
        BLV_minus1(C==1)=E;
        BLV_minus1(C==2)=0;
        BLV_minus1(BLV_minus1<0)=0;
        % assume SS make this the ratio of densities
    else
        BLV_minus1 = BLV;
    end
    
    S_minus1 = B+H;
    
    Hbeforesoilprod = H;
    
	% exponential soil production
    [H, B, BLV, X, Ns_soil_prod] = SoilProd(H,B,X,C,nX,alpha,P,rhos,rhor,xr,dt,M2,M,E,Ns,Nzb,host_min);
    
    BLV(BLV<0)=0;
    BLV_minus1(BLV_minus1<0)=0;    
    
    zr_rat = X(:,:,zr_min)./xr(zr_min); % zirconium concentration ratio
    zr_rat_inv = xr(zr_min)./X(:,:,zr_min);
    
    % Alternate calculation stemming from the solution of the mass balance
    % PDE for Ns using both Ns and Nzb -- highly accurate but not practical for field workers as
	% two concentrations needed
    Dalt_prod = (cosmo_prod_spal.*L1.*(1 - exp(-rhos.*H/L1))./(rhos.*H)) + (cosmo_prod_neg_1.*L2.*(1 - exp(-rhos.*H/L2))./(rhos.*H)) + (cosmo_prod_fast.*L4.*(1 - exp(-rhos.*H/L4))./(rhos.*H));
    Dalt_m = (Dalt_prod - Ns*lambda_10Be).*(rhos.*H).*zr_rat./(Ns - Nzb);
    Dalt = Dalt_m./rhos;
    
	% this routine gives the mean BLV for filling in the bottom of the TCN profile comment out
	% if using cosmo_noprof()
    if n==1
        BLV_cube=BLV_cube; % time-cube of blv rates
        BLV_mat=reshape(mean(BLV_cube),[K,J]); % mean of cube shaped to landscape grid
    else
        BLV_cube=circshift(BLV_cube,1); % circularly shift values by 1 
        BLV_cube(1,:,:)=BLV;
        BLV_mat=reshape(mean(BLV_cube),[K,J]);
    end   
    
    S=B+H; % surface for nuclide production
    
    % Comment or uncomment, depending on whether you want variable
    % production rates
    cosmo_prod_spal =  ppval(P_spal_pchip,S)*1000; % spal production for 10Be (atoms kg-1 yr-1) at surface (Balco et al., 2008) (p) 
    cosmo_prod_neg_1 = ppval(P_neg_pchip,S)*1000; % first neg moun production rate (modeled as two exponentials) (p)
    cosmo_prod_fast = ppval(P_fast_pchip,S)*1000;
     
	% comment out of the these
	% semi-Lagrangian advection with bedrock profile TCN profile tracking (slow) 
    % [Ns, Nzb, N_prof, N_prof2, Dinf_m, Dinf_m_lal, Dinf_m_zb] = cosmo_prof(BLV, BLV_minus1, BLV_mat, H, N_prof,...
    % N_prof2, Ns, Nzb, C, rhor, rhos, depth, depth_resolution, cosmo_prod_spal, cosmo_prod_neg_1,...
    % cosmo_prod_fast, lambda_10Be, L1, L2, L4, dt, t, K, J, xr, X, P, alpha,...
    % M2, M, host_min, zr_min);
    
	% No profile tracking -- pseudo-profile based on current Nzb (fast)
	[Ns, Nzb, Nzb_spal, Nzb_neg, Nzb_fast, Dinf_m, Dinf_m_lal, Dinf_m_zb] = cosmo_noprof(H, Ns, Nzb, Nzb_spal, Nzb_neg, Nzb_fast, C, rhor, rhos, cosmo_prod_spal, cosmo_prod_neg_1,...
    cosmo_prod_fast, lambda_10Be, L1, L2, L4, dt, t, K, J, xr, X, P, alpha, M2, M, host_min, zr_min)
	
	% No Ns in channels
	Ns(C==1)=0;
    
    % Mineral dissolution based on concentration and dissolution rate constants
    [H, X] = Weather(H,X,C,nX,kx,Ax,sx,wx,rhos,xr,dt,M2,M);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%% Boundary conditions/perturbations %%%%%%%%%%%%%%

    % Stream incision/relative uplift
    
    % For lake infilling and stream incision cessation, no uplift is needed
    % anywhere
    
    S = B + H;
    
    B(C~=1) = B(C~=1) + E*dt;
    
	% Rock uplift/relative base-level lowering at boundaries 
	% The more complex routine below changes channel elevation to match channel-side denudation rates
	% so the channel does not rise above the hillslopes in the climatic scenario -- this adjustment is slight
    B(M==1) = B(M==1) - E*dt;
    
    for i=1:K-1
      for j=1:J-1
          if C(i,j)==1 & B(i,j)>=S(i+1,j) & C(i+1,j)~=1
                  B(i,j)=S(i+1,j);
          elseif C(i,j)==1 & B(i,j)>=S(i,j+1) & C(i,j+1)~=1
                  B(i,j)=S(i,j+1);
          end
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%% Surface Evolution %%%%%%%%%%%%%%%%%%%%%%%
   
    Hbeforeerosion = H; % use to infer erosion rate grid
    Xbeforeerosion = X;
	
    [H, B, X, Ns] = Erode_X_NS(H,B,X,Ns,Nzb,C,nX,K,J,dt,dx,dy,S_c,K3,xr,M,M2,rhos,host_min);
    Haftererosion = H;
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%%%%%%%%%%%%%%%%%%% Store previous timestep %%%%%%%%%%%%%%%%%%%%

    Hnminus1 = Htemp;
    Bnminus1 = Btemp;
    Xnminus1 = Xtemp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    %%%%%%%%%% Physical and chemical erosion rates %%%%%%%%%%%%%%%%%
    
    % Compute grids of physical (R) and chemical (W) erosion rates (m/yr).  
    % Note that these are thinning rates of the soil, so both would have to
    % be multiplied by rhos to obtain the mass flux out of the soil.

    % Compute W as sum of chemical erosion fluxes from all mineral phases.
    temp = zeros(size(X));
    for i = 1:length(kx)
        temp(:,:,i) = kx(i)*Ax(i)*X(:,:,i) - sx(i)*wx(i)/rhos;
    end
    W = H .* sum(temp,3);  % Sum over all mineral phases (m/yr)
    clear temp
    
    % Physical erosion and denudation (instantaneous) 
    R = (Hbeforeerosion - Haftererosion)/dt;  % m/yr
    R_m = R*rhos; % erosion in mass form
    
    if n == 1
        Dact_old = output.Dact_mperyr(:,:,1); % denudation rate (n-1) [m yr-1]
    else
        Dact_old = Dact;
    end
    
    % this will actually be Dact at n-1 when a proper steady state is
    % avail.
    Dact = W + R; % bedrock lowering denudation rates [m yr-1]
    Dact(C==1) = 0;
    Dact(C==2) = 0; % boundaries
    Dact(M==1) = 0;
    R(M==1)=0;
    W(M==1)=0;
    
    
    S=B+H;
    
 
    
    Dinf_m(C==1) = 0;
    Dinf_m(C==2) = 0;
    Dinf = Dinf_m./rhos;
    Dinf_lal = Dinf_m_lal./rhos;
    Dinf_zb = Dinf_m_zb./rhos;
    D_diff = (Dact - Dinf);
    CDF2 = 1 - xr(zr_min)./Xbeforeerosion(:,:,zr_min);
    Winf = Dinf.*CDF2;
    Einf = Dinf - Winf;
    Winf_alt = Dalt.*CDF2;
    Einf_alt = Dalt - Winf_alt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if do_output  % only if we're saving results
        if rem(n,save_interval)==0
            lastsave = n/save_interval+1;  % Start at second save timestep.
            output.topo_m(:,:,lastsave) = B + H;  % topography (m)
            output.H_m(:,:,lastsave) = H;  % soil thickness (m)
            output.X(:,:,:,lastsave) = X;  % soil mineral abundances ()
            output.W_mperyr(:,:,lastsave) = W;  % chemical erosion rate (m/yr)
            output.R_mperyr(:,:,lastsave) = R;  % physical erosion rate (m/yr)
            output.Dact_mperyr(:,:,lastsave) = Dact;
            output.Dinf_mperyr(:,:,lastsave) = Dinf; % inferred denudation rate (m/yr)
            output.Dinf_lal_mperyr(:,:,lastsave) = Dinf_lal;
            output.Dinf_zb_mperyr(:,:,lastsave) = Dinf_zb;
            output.Dalt_mperyr(:,:,lastsave) = Dalt;
            output.Winf_mperyr(:,:,lastsave) = Winf;
            output.Einf_mperyr(:,:,lastsave) = Einf;
            output.Einf_alt_mperyr(:,:,lastsave) = Einf_alt;
            output.Winf_alt_mperyr(:,:,lastsave) = Winf_alt;
            output.Ddiff_mperyr(:,:,lastsave) = D_diff;
            output.Ns_atg(:,:,lastsave) = Ns/1000;
            output.Nzb_atg(:,:,lastsave) = Nzb/1000;
            %output.Ns_lake(lastsave) = Ns_lake;
            output.t_yr(lastsave) = n*dt;  % time since onset of simulation (yr)
            output.t_yr(lastsave) = t;  % time since onset of simulation (yr)
        end
    end
    

	% This for tailing output at certain grid cell	-- can actually slow down runs when not tracking TCN
	% profiles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('---------------------------------------------\n')
    fprintf('timestep: %d\n', n) 
    fprintf('inferred denudation rate at (58,74): %s\n', Dinf(58,74))
    fprintf('inferred denudation rate [Lal] at (58,74): %s\n', Dinf_lal(58,74))
    fprintf('inferred denudation rate [Bedrock] at (58,74): %s\n', Dinf_zb(58,74))
    fprintf('actual denudation rate at (58,74): %s\n', Dact(58,74))
    fprintf('B at (58,74): %s\n', B(58,74))
    
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END MAIN ITERATION LOOP % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute CDF

CDF = zeros(size(output.H_m));
zirconindex = find(string(p.xname)=='zircon');
for i = 1:size(output.X,4)
    CDF(:,:,i) = 1 - p.xr(zirconindex)./output.X(:,:,zirconindex,i);
end


% Set edges to NaNs for certain variables.
CDF(C==1) = NaN;
output.H_m(C==1) = NaN;
output.R_mperyr(C==1) = NaN;
output.W_mperyr(C==1) = NaN;
output.X(C==1) = NaN;
output.Dact(C==1) = NaN;
output.Dinf(C==1) = NaN;


%% Save inputs and outputs.

if do_output % only if we're saving results
    outputname = ['output_' run_name '_' clocktime2text '.mat'];
    save(outputname, 'inputs', 'p', 'output');
end


