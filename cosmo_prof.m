% Cosmo operations where getting pretty lengthy so I am making a separate
% m-file....This is not part of the previous operator-splitting scheme, but
% its own forward scheme with various sub-schemes(?). Still needs chemical
% erosion and decay. 


function [Ns, Nzb, N_prof, N_prof2, Dinf_m, Dinf_m_lal, Dinf_m_zb] = cosmo_prof(BLV, BLV_minus1, BLV_mat, H, N_prof,...
    N_prof2, Ns, Nzb,C, rhor, rhos, depth, depth_resolution, cosmo_prod_spal, cosmo_prod_neg_1,...
    cosmo_prod_fast, lambda_10Be, L1, L2, L4, dt, t, K, J, xr, X, P, alpha, M2, M, host_min, zr_min)

% This will be integrated into an .m file for cosmo production


H2 = H; %Hbeforeerosion; % depth of soil above bedrock surface (Hbeforeerosion)
depth_ints = linspace(0, depth, depth_resolution); % (p) (i.e., 200 equally space pts below surface)
depth_ints = depth_ints'.*ones(length(depth_ints),K,J);
depth_ints = permute(depth_ints, [2,3,1]);
depth2 = 2; % lower depth in meters
depth_resolution2 = 100; % coarser resolution for deeper production (test)
depth_ints2 = linspace(depth, depth2, depth_resolution2);
depth_ints2 = depth_ints2'.*ones(length(depth_ints2),K,J);
depth_ints2 = permute(depth_ints2, [2,3,1]);
H_mat = repmat(H2,1,1,depth_resolution); % add another dimension to H with same value at every belowground point
H_mat2 = repmat(H2,1,1,depth_resolution2); 
depth_ints_r = depth_ints+H_mat; % working
depth_ints_r2 = depth_ints2+H_mat2;
md_soil = H_mat*rhos; % mass-depth (kg m-2) of soil (solved each step) -- which H? Ns was solved before Erode() should it now fall after SoilProd()?
md_soil2 = H_mat2*rhos;
md_rock = depth_ints_r*rhor - md_soil; % mass-depth of rock down profile (solved each step)
md_rock2 = depth_ints_r2*rhor - md_soil2;
V = BLV_minus1; % bedrock lowering velocity at n-1
Vnext = BLV; % bedrock lowering velocity at n [m yr-1]
clear depth_ints depth_ints2 depth_ints_r depth_ints_r2 H_mat H_mat2

% Granger and Riebe eqns for production at depth (2014)
% (solved every step to find source terms [points won't 'know' whey they
% are within the rock column but source term will])


% Shielding parameter could be inserted here 
% Fits obtained from a combination of Balco's and Marrero's codes
% 

%Ps_S=fit_LSD_Ps.coeff(1).*S.^4 + fit_LSD_Ps.coeff(2).*S.^3 + fit_LSD_Ps.coeff(3).*S.^2 + ...
%     fit_LSD_Ps.coeff(4).*S + fit_LSD_Ps.coeff(5);
 
%Pneg_S=fit_LSD_neg.coeff(1).*S.^4 + fit_LSD_neg.coeff(2).*S.^3 + fit_LSD_neg.coeff(3).*S.^2 + ...
%     fit_LSD_neg.coeff(4).*S + fit_LSD_neg.coeff(5);
 
%Pfast_S=fit_LSD_fast.coeff(1).*S.^4 + fit_LSD_fast.coeff(2).*S.^3 + fit_LSD_fast.coeff(3).*S.^2 + ...
%     fit_LSD_fast.coeff(4).*S + fit_LSD_fast.coeff(5);

%Ps = Ps_S.*exp(-md_rock./L1);



Ps = cosmo_prod_spal.*exp(-md_rock./L1);
Pn_1 = cosmo_prod_neg_1.*exp(-md_rock./L2);
%Pn_2 = cosmo_prod_neg_2*exp(-md_rock/L3);
Pf = cosmo_prod_fast.*exp(-md_rock./L4);
Ptot = (Ps + Pn_1 + Pf);
clear Ps Pn_1 Pf

Ps_2 = cosmo_prod_spal.*exp(-md_rock2./L1);
Pn_1_2 = cosmo_prod_neg_1.*exp(-md_rock2./L2);
%Pn_2_2 = cosmo_prod_neg_2*exp(-md_rock2/L3);
Pf_2 = cosmo_prod_fast.*exp(-md_rock2./L4);
Ptot_2 = (Ps_2 + Pn_1_2 + Pf_2);
clear Ps_2 Pn_1_2 Pn_2_2 Pf_2

%Ptot(M_cube==1)=NaN;
%Ptot_2(M_cube2==1)=NaN;


% Semi-Lagrangian advection scheme to keep track of deep cosmogenic
% profiles based on code by Blackburn et al. (2018) for crustal temp profile. 
% Muogenic production approximated by a series of exponentials (Granger and Muzikar, 2001).
% Ptot is total production of nuclides at different depths
% D_old is the n-1 denudation rate (physical + chemical)
% D is the n denudation rate
%erosion rate (m/yr) (supplied to scheme Hbeforeerosion_old - Hbeforeerosion; gotta get this right)

% Lower column(s) must solved first then handed off to upper column

J3 = length(Ptot_2); % eventually fix all these variable names to reflect what is going on
jvec2 = flip(1:J3)';
jvec2 = repmat(jvec2, 1, K, J);
jvec2 = permute(jvec2, [2, 3, 1]);
S2 = Ptot_2; % production at depth intervals --> source terms
dz2 = depth/depth_resolution; % space between points on vertical grid (m) depth/depth_resolution
dt = dt; % p.dt
% Make V, Vnext, S vectors of size(U) if they aren't already
V2 = V.*ones(K,J,depth_resolution2); 
Vnext2 = Vnext.*ones(K,J,depth_resolution2);

% number of iterations to do when finding departure point
numit = 2; 

% Interpolate velocities half a time step in the future
Vhalf = 0.5*(V2 + Vnext2); % has to do with the flavor of semi-Langragian scheme 

% Set initial guess velocity to velocity at time n+1
Vtemp2 = Vnext2;

[X3,Y2,Z2] = meshgrid(1:K, 1:J,1:depth_resolution2);
A2 = [2 1 3];
X3 = permute(X3,A2);
Y2 = permute(Y2,A2);
Z2 = permute(Z2,A2);
V2 = permute(V2,A2);

for k=1:numit

    % get fractional j values half a time step in the future
    %jmid = jvec+0.5*dt*Vtemp/dz; 
    Z2a = Z2+0.5*dt*Vtemp2/dz2;
    
    % interpolate linearly to get velocities half a time step in the future
    %Vtemp = interp1(jvec,Vhalf,jmid,'linear'); 
    F1_2 = griddedInterpolant(X3,Y2,Z2,Vhalf,'linear','none');
    Vtemp2 = F1_2(X3,Y2,Z2a);
end

% get departure points at present time
%jdep = jvec+dt*Vtemp/dz; 
Z2b = Z2+dt*Vtemp2/dz2;
% Interpolate to find source term value at departure points
%Sdep = interp1(jvec,S,jdep,'PCHIP'); 
%S2(M_cube2==1)=0.001;
F2_2 = griddedInterpolant(X3,Y2,Z2,S2,'cubic','none');
%Sdep2 = nan(size(S2));
%Sdep2(M_cube2==0) = F2_2(X3(M_cube2==0),Y2(M_cube2==0),Z2b(M_cube2==0));
Sdep2 = F2_2(X3,Y2,Z2b);

% 'Get departure elevations by interpolation, and add source term
% Use trapezoidal rule, a la Spiegelman, to estimate the integral of the 
% source term along the solution characteristic (i.e., average source term 
% values at departure and future points. We do not account for 
% time-dependence of source term.' --Blackburn

%N_prof2(M_cube==1)=0.001;
F3_2 = griddedInterpolant(X3,Y2,Z2,N_prof2, 'cubic','none');
%T=interp1(jvec,N_prof,jdep,'makima')+ dt*0.5*(Sdep+S); % S * dt  = K/s * s = K
%N_prof_new2 = nan(size(N_prof2));
%N_prof_new2(M_cube2==0) = F3_2(X3(M_cube2==0),Y2(M_cube2==0),Z2b(M_cube2==0)) + dt*0.5*(Sdep2(M_cube2==0)+S2(M_cube2==0)) ;
N_prof_new2 = F3_2(X3,Y2,Z2b);
N_prof_new2 = N_prof_new2 - lambda_10Be*N_prof_new2*dt;

% Reduce profile concentrations due to decay and set the interpolated
% profile to N_prof
N_prof2 = N_prof_new2 + dt*0.5*(Sdep2+S2);

% Lowermost point on vertical grids needs filled
Nspal_last = cosmo_prod_spal.*exp(-md_rock2./L1)./(lambda_10Be + (BLV_mat*rhor)/L1); %BLV_mat
Nn_1_last = cosmo_prod_neg_1.*exp(-md_rock2./L2)./(lambda_10Be + (BLV_mat*rhor)/L2);
%Nn_2_last = cosmo_prod_neg_2.*exp(-md_rock2/L3)./(lambda_10Be + (BLV_mat*rhor)/L3);
Nf_last = cosmo_prod_fast.*exp(-md_rock2./L4)./(lambda_10Be + (BLV_mat*rhor)/L4);
Ntot_last = (Nspal_last + Nn_1_last + Nf_last);
N_prof2(:,:,end)=Ntot_last(:,:,end);  %testing

J2 = length(Ptot); % eventually fix all these variable names to reflect what is going on
jvec = flip(1:J2)';
jvec = repmat(jvec, 1, K, J);
jvec = permute(jvec, [2, 3, 1]);
S = Ptot; % production at depth intervals --> source terms
dz = depth/depth_resolution; % space between points on vertical grid (m) depth/depth_resolution
dt = dt; % p.dt
% Make V, Vnext, S vectors of size(U) if they aren't already
V = V.*ones(K,J,depth_resolution); 
Vnext = Vnext.*ones(K,J,depth_resolution);

% number of iterations to do when finding departure point
numit = 2; 

% Interpolate velocities half a time step in the future
Vhalf = 0.5*(V + Vnext); % has to do with the flavor of semi-Langragian scheme 

% Set initial guess velocity to velocity at time n+1
Vtemp = Vnext;

[X2,Y,Z] = meshgrid(1:K, 1:J,1:depth_resolution);
A = [2 1 3];
X2 = permute(X2,A);
Y = permute(Y,A);
Z = permute(Z,A);
V = permute(V,A);

for k=1:numit

    % get fractional j values half a time step in the future
    %jmid = jvec+0.5*dt*Vtemp/dz; 
    Z2 = Z+0.5*dt*Vtemp/dz;
    % interpolate linearly to get velocities half a time step in the future
    %Vtemp = interp1(jvec,Vhalf,jmid,'linear'); 
    F1 = griddedInterpolant(X2,Y,Z,Vhalf,'linear','none');
    Vtemp = F1(X2,Y,Z2);
end

% get departure points at present time
%jdep = jvec+dt*Vtemp/dz; 
Z3 = Z+dt*Vtemp/dz;
% Interpolate to find source term value at departure points
%Sdep = interp1(jvec,S,jdep,'PCHIP'); 
%S(M_cube==1)=0.001;
F2 = griddedInterpolant(X2,Y,Z,S,'cubic','none');
%Sdep = nan(size(S));
%Sdep(M_cube==0) = F2(X2(M_cube==0),Y(M_cube==0),Z3(M_cube==0)); 
Sdep = F2(X2,Y,Z3);
% 
% % Get departure elevations by interpolation, and add source term
% % Use trapezoidal rule, a la Spiegelman, to estimate the integral of the 
% % source term along the solution characteristic (i.e., average source term 
% % values at departure and future points. We do not account for 
% % time-dependence of source term. --Blackburn
% 
%N_prof(M_cube==1)=0.001;
F3 = griddedInterpolant(X2,Y,Z,N_prof, 'cubic','none');
%T=interp1(jvec,N_prof,jdep,'makima')+ dt*0.5*(Sdep+S); % S * dt  = K/s * s = K
%N_prof_new = nan(size(N_prof));
%N_prof_new(M_cube==0) = F3(X2(M_cube==0),Y(M_cube==0),Z3(M_cube==0)) + dt*0.5*(Sdep(M_cube==0)+S(M_cube==0)) ;
N_prof_new = F3(X2,Y,Z3);
% Reduce profile concentrations due to decay and set the interpolated
% profile to N_prof
N_prof_new = N_prof_new - lambda_10Be*N_prof_new*dt ;
N_prof = N_prof_new + dt*0.5*(Sdep+S);

N_prof(:,:,end)=N_prof2(:,:,1); 

clear N_prof_new N_prof_new2
% 
% 
if sum(isnan(Nzb(:)))>0
    fprintf('cosmo.m: Nzb NaN(s) detected before timestep \n')
end
if sum(isnan(Ns(:)))>0
    fprintf('cosmo.m: Ns NaN(s) detected before timestep \n')
end

if sum(isnan(N_prof(:)))>0
    fprintf('cosmo.m: N_prof NaN(s) detected before timestep \n')
end

if sum(isnan(N_prof2(:)))>0
    fprintf('cosmo.m: N_prof2 NaN(s) detected before timestep \n')
end
% 
% 
% 
% 
% %Nzb = Nzb - lambda_10Be*Nzb*dt;
% 
Ns_spal = ((cosmo_prod_spal.*L1.*(1-exp(-rhos.*H./L1))));
Ns_neg_1 = ((cosmo_prod_neg_1.*L2.*(1-exp(-rhos.*H./L2)))); 
% Ns_neg_2 = ((cosmo_prod_neg_2*L3.*(1-exp(-rhos.*H/L3))));
Ns_fast = ((cosmo_prod_fast.*L4.*(1-exp(-rhos.*H./L4))));
Ns_prod_tot = (Ns_spal + Ns_neg_1 + Ns_fast)./(rhos.*H); % If using two expoentials for neg muons add here
% 
% 
% 
% 
% % The change in Ns due to hillslope transport is in Erode.m
% 
% %fprintf('cosmo - soil 10Be concentration at (200,200): %s\n', Ns(200,200))
% %fprintf('cosmo - H at (200,200): %s\n', H(200,200))
% 
Ns_1 = Ns + Ns_prod_tot*dt - ((xr(host_min).*Ns*P.*exp(-alpha.*H))./(rhos.*X(:,:,host_min).*H))*dt + ((xr(host_min).*Nzb*P.*exp(-alpha.*H))./(rhos.*X(:,:,host_min).*H))*dt - lambda_10Be*Ns*dt;
% %Ns_1 = Ns + Ns_prod_tot*dt - Ns_soil_prod - lambda_10Be*Ns*dt;
% 
Ns(M==0)=Ns_1(M==0); % Everything outside of catchment remains the same

Nzb = N_prof(:,:,1);
% 
% Zero out soil in channels, boundaries, and lake cells 

Ns(M2==1)=0;
Ns(C==1) = 0;
Ns(C==2) = 0;
Ns(M==1)=0;

% 
% % Zero out bedrock in channels, boundaries, and lake cells
Nzb(M2==1)=0;
Nzb(C==1)=0;
Nzb(C==2)=0;
% 
% 
% 
% 
% % Ns cannot be zero -- very high fluxes + 
% 
% %Ns(Ns<Nzb)=Nzb(Ns<Nzb);
% 
% 
% if sum(isnan(Ns(:)))>0
%     fprintf('cosmo.m: Ns NaN(s) detected after timestep \n')
% end
% 
% %fprintf('cosmo - soil 10Be concentration at (200,118): %s\n', Ns(200,118))
% 
zr_rat = X(:,:,zr_min)./xr(zr_min);

Dinf_spal = (cosmo_prod_spal.*L1./Ns).*(zr_rat.*(1-exp(-rhos.*H./L1)) + exp(-rhos.*H./L1)); 
Dinf_neg_1 = (cosmo_prod_neg_1.*L2./Ns).*(zr_rat.*(1-exp(-rhos.*H./L2)) + exp(-rhos.*H./L2)); 

%Dinf_neg_2 = (cosmo_prod_neg_2*L3./Ns).*((X(:,:,2)./xr(2).*(1-exp(-rhos.*H2/L3))+exp(-rhos.*H2/L3))) - lambda_10Be*rhor*H2; 
Dinf_fast = (cosmo_prod_fast.*L4./Ns).*(zr_rat.*(1-exp(-rhos.*H./L4)) + exp(-rhos.*H./L4));
Dinf_m = Dinf_spal + Dinf_neg_1 + Dinf_fast;


Dinf_spal_lal = (cosmo_prod_spal.*L1./Ns);
Dinf_neg_1_lal = (cosmo_prod_neg_1.*L2./Ns);

%Dinf_neg_2 = (cosmo_prod_neg_2*L3./Ns).*((X(:,:,2)./xr(2).*(1-exp(-rhos.*H2/L3))+exp(-rhos.*H2/L3))) - lambda_10Be*rhor*H2; 
Dinf_fast_lal = (cosmo_prod_fast.*L4./Ns);
Dinf_m_lal = Dinf_spal_lal + Dinf_neg_1_lal + Dinf_fast_lal;

Dinf_spal_zb = (cosmo_prod_spal.*L1./Nzb).*exp(-rhos.*H./L1);
Dinf_neg_zb = (cosmo_prod_neg_1.*L2./Nzb).*exp(-rhos.*H./L2);
Dinf_fast_zb = (cosmo_prod_fast.*L4./Nzb).*exp(-rhos.*H./L4);
Dinf_m_zb = Dinf_spal_zb + Dinf_neg_zb + Dinf_fast_zb;

Dinf_m(C==1)=0;
Dinf_m(C==2)=0;
Dinf_m(M2==1)=0;
Dinf_m(M==1)=0;

Dinf_m_lal(C==1)=0;
Dinf_m_lal(C==2)=0;
Dinf_m_lal(M2==1)=0;
Dinf_m_lal(M==1)=0;

Dinf_m_zb(C==1)=0;
Dinf_m_zb(C==2)=0;
Dinf_m_zb(M2==1)=0;
Dinf_m_zb(M==1)=0;

% Short-circuit to eliminate function 

% Ns = Ns;
% Nzb = Nzb;
% N_prof = N_prof;
% N_prof2 = N_prof2;
% Dinf_m = zeros(size(Ns));

end