% Params
p.dx = 10; % grid spacing x-direction [m] 
p.dy = 10; %   "    "     y-direction [m]
p.rhor = 2650; % bedrock density [kg m-3] - constant
p.rhos = p.rhor/2; % soil density [kg m-3] - constant
p.E = ; % bedrock uplift rate [m yr-1]
p.D = 0.003; % nonlinear hillslope diffusion coefficient [m2 yr-1]
p.P = 0.7; % soil production coefficient/max production rate [kg m-2 yr-1]
p.alpha = 3; 
p.nX = 2; % # minerals in bedrock -- can be any number; just would need extra parameterization
p.kx = [5e-06; 0]; % mineral reactivity [mol m-2 yr-1]
p.Ax = [117 4]; % mineral specific surface area [m2 mol-1]
p.sx = [0 0]; % clay production rate [mol m-3 yr-1]
p.wx = [0.26544 0.18331]; % mineral molar mass [kg mol^-1]
p.xr = [0.5 0.5]; % bedrock mineral concentrations [kg kg-1]
p.xname = {'plagioclase' 'zircon'}; % mineral names
p.host_min = 2; % cosmogenic host mineral
p.zr_min = 2; % inert mineral for denudation rate paritioning
p.K = 150; % domain size x-direction [# cells]
p.K = 150; % domain size y-direction [# cells]

% C - Channel network (channel network == 1; hillslope == 0) | M - boundaries (boundary == 1; non-boundary == 0)
% Use your prefered drainage area algorithm with a channel threshold area to extract the grid
% Make sure there are no cells where C == 1 & M == 1 or you get double incision

% If static points are needed
M2 = zeros(size[K,J]); 

% Generate X (grid of mineral concentrations)
x1 = [0.3 0.7]; % initial guess for steady-state soil concentrations

for i=p.nX
    Xa=zeros(size([K,J]))+x1(i);
    X(:,:,i)=Xa;
end

for i=1:K
    for j=1:J
        if C(i,j)==1
            X(i,j,:)=p.xr;
        end
    end
end

for i=1:K
    for j=1:J
        if M(i,j)==1
            X(i,j,:)=p.xr;
        end
    end
end

H = zeros(size([K,J]))+0.5; % initial guess at steady-state soil thickness 

% B is bedrock elevation. You generate this with another landscape evolution model or programatically
S = B + H; % surface elevation [m] 

% production pchips from var_prod.m (requires cronuscalc to scale) or load pchips.mat
cosmo_prod_spal =  ppval(P_spal_pchip,S)*1000; % spal production for 10Be (atoms kg-1 yr-1) at surface (Balco et al., 2008) (p) 
cosmo_prod_neg_1 = ppval(P_neg_pchip,S)*1000; % first neg moun production rate (modeled as two exponentials) (p)
cosmo_prod_fast = ppval(P_fast_pchip,S)*1000; % fast muon production rate (p) 


% Mean attenuation lengths determined by Balco's muon codes -- Loaded from pchips.mat or cronuscalc or defined here
% L1 = 1522.8; % spal atten. (kg m-2) (many); ranges between 1500-1900 (site-calibrated)
% L2 = 13675; % neg
% L4 = 35144; % fast muon 

host_min = p.host_min; % host mineral #; refer to p.xname
zr_min = p.zr_min; % zircon mineral #; for Dinf & CDF


t_half_10Be = 1.387e06; % Chemeleff et al., 2010
lambda_10Be = log(2)/t_half_10Be; % 10Be decay constant

p.lambda_10Be = lambda_10be;

K = p.K;
J = p.J;

% For semi-Langrangian advection scheme
depth = 1; % depth of profile #1 [m]
depth_resolution = 100; % # points in profile #1 - 1/100 = 100 pts m-1
depth_ints = linspace(0, depth, depth_resolution);
depth_ints = depth_ints'.*ones(length(depth_ints),K,J);
 
depth2 = 2; % depth of profile #2 [m] - The resolutions of profile #2 can be different (i.e., deeper --> lower res)
depth_resolution2 = 100;
depth_ints2 = linspace(depth, depth2, depth_resolution2);
depth_ints2 = depth_ints2'.*ones(length(depth_ints2),K,J);

% Create profiles of mass-depth
depth_ints = permute(depth_ints, [2,3,1]);
H_mat = repmat(H,1,1,depth_resolution); % add another dimension to H with same value at every belowground point
depth_ints_r = depth_ints+H_mat; % working
md_soil = H_mat*rhos; % mass-depth (kg m-2) of soil (solved each step) -- which H? Ns was solved before Erode() should it now fall after SoilProd()?
md_rock = depth_ints_r*rhor - md_soil; 

depth_ints2 = permute(depth_ints2, [2,3,1]);
H_mat2 = repmat(H,1,1,depth_resolution2); % add another dimension to H with same value at every belowground point
depth_ints_r2 = depth_ints2+H_mat2; % working
md_soil2 = H_mat2*rhos; % mass-depth (kg m-2) of soil (solved each step) -- which H? Ns was solved before Erode() should it now fall after SoilProd()?
md_rock2 = depth_ints_r2*rhor - md_soil2; 

clear depth_ints depth_ints2 depth_ints_r depth_ints_r2 H_mat H_mat2 md_soil md_soil2

% This matrix is for averaging the bedrock-lowering rate in transient
% scenarios (boundary filling routine) -- this could be better

init_BLV = p.E; % bedrock-lowering velocity (m yr-1)
BLV = init_BLV*(ones(K,J));
BLV_minus1 = BLV;
BLV_cube=init_BLV*ones(200,K,J); % (blv @ time, X, Y) the cube will be used to create 2D matrix of 
% time-averaged bedrock-lowering values

% soil denudation has to equal total denudation (i.e., no saprolite
% weathering -- Granger and Riebe (2014) eqn. 19 for Nsap but should work
% for bedrock profile as well (solved before main loop to instantiate
% concentration profile) -- not needed if steady state profiles exists

initial_Dmass = init_BLV*p.rhor; % denudation rate as mass [kg m-2 yr-1] 

Nsp = cosmo_prod_spal.*exp(-md_rock/L1)/(lambda_10Be + initial_Dmass/L1);
Nn_1 = cosmo_prod_neg_1.*exp(-md_rock/L2)/(lambda_10Be + initial_Dmass/L2);
Nf = cosmo_prod_fast.*exp(-md_rock/L4)/(lambda_10Be + initial_Dmass/L4);

N_prof = (Nsp + Nn_1 + Nf);

% If not using the semi-Lagrangian advection scheme (no TCN profile tracking) we need the fraction
% of Nzb from each production pathway (31337 mode)

Nzb_spal = Nsp(:,:,1); % spallation
Nzb_neg = Nn_1(:,:,1); % negative muon capture
Nzb_fast = Nf(:,:,1); % fast muons

clear Nsp Nn_1 Nf

Nsp_2 = cosmo_prod_spal.*exp(-md_rock2/L1)/(lambda_10Be + initial_Dmass/L1);
Nn_1_2 = cosmo_prod_neg_1.*exp(-md_rock2/L2)/(lambda_10Be + initial_Dmass/L2);
Nf_2 = cosmo_prod_fast.*exp(-md_rock2/L4)/(lambda_10Be + initial_Dmass/L4);
N_prof2 = (Nsp_2 + Nn_1_2 + Nf_2);

clear Nsp_2 Nn_1_2 Nf_2

% Bedrock conc. is the top profile point
Nzb = N_prof(:,:,1);



% solve for steady-state concentration in soil if none exists

Ns_spal = cosmo_prod_spal./(lambda_10Be + initial_Dmass./L1).*((X(:,:,host_min)./xr(host_min).*(1-exp(-rhos.*H./L1))+exp(-rhos.*H./L1)));
Ns_neg_1 = cosmo_prod_neg_1./(lambda_10Be + initial_Dmass./L2).*((X(:,:,host_min)./xr(host_min).*(1-exp(-rhos.*H./L2))+exp(-rhos.*H./L2)));
Ns_f = cosmo_prod_fast./(lambda_10Be + initial_Dmass./L4).*((X(:,:,host_min)./xr(host_min).*(1-exp(-rhos.*H./L4))+exp(-rhos.*H./L4)));
Ns = Ns_spal + Ns_neg_1 + Ns_f;
clear Ns_spal Ns_neg_1 Ns_f

% No concentrations in static channels or outside of the domain
% Ns(C==1)=0;
% Ns(M2==1)=0;



%% initial inferred denudation rate and associated inferred physical and chemical erosion rates

Dinf_spal = (cosmo_prod_spal.*L1./Ns).*((X(:,:,zr_min)./xr(zr_min).*(1-exp(-rhos.*H./L1))+exp(-rhos.*H./L1))); 
Dinf_neg_1 = (cosmo_prod_neg_1.*L2./Ns).*((X(:,:,zr_min)./xr(zr_min).*(1-exp(-rhos.*H./L2))+exp(-rhos.*H./L2))); 
%Dinf_neg_2 = (cosmo_prod_neg_2*L3./Ns).*((X(:,:,zr_min)./xr(zr_min).*(1-exp(-rhos.*H/L3))+exp(-rhos.*H/L3))) - lambda_10Be*rhor.*H; 
Dinf_fast = (cosmo_prod_fast.*L4./Ns).*((X(:,:,zr_min)./xr(zr_min).*(1-exp(-rhos.*H./L4))+exp(-rhos.*H./L4)));
Dinf_m = Dinf_spal + Dinf_neg_1 + Dinf_fast;
Dinf = Dinf_m/rhos; 
Dinf(C==1)=0;
Dinf(M==1)=0;
output.Dinf_mperyr(:,:,1) = Dinf;
CDF_start = 1 - xr(zr_min)./X(:,:,zr_min);
Winf = Dinf.*CDF_start;
Einf = Dinf - Winf;
output.Winf_mperyr(:,:,1) = Winf;
output.Einf_mperyr(:,:,1) = Einf;