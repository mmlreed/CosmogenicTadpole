function [H, B, BLV, X, Ns_soil_prod] = SoilProd(H,B,X,C,nX,alpha,P,rhos,rhor,xr,dt,M2,M,E,Ns,Nzb,host_min)
%
% Update soil thickness H, bedrock elevation B, mineral abundances X, and
% cosmogenic radionuclide concentrations in bedrock and soil according to 
% the perturbations to these quantities by soil production.
% This is part of a splitting method, whereby this function computes a
% portion of the temporal changes in H, B, and X, and other functions
% compute the other portions of those changes.


%% Update soil thickness H
% Use exact solution for change in soil thickness over an interval dt.
% This solution can be derived via separation of variables.
H(C==1) = 1e-10;
H(C==2) = 1e-10;
H(M2==1)=1e-10;
H(M==1)=1e-10;
Hnplus1 = (1/alpha)*log(alpha*P*dt/rhos + exp(alpha*H));

if sum(isnan(H(:)))>0
    fprintf('SoilProd.m: NaN(s) detected in H after soil production\n')
end

%fprintf('soil_prod - soil 10Be concentration at (200,200): %s\n', Ns(200,200))
%fprintf('soil_prod - H at (200,200): %s\n', H(200,200))

% This is to limit extremely thin soils
Hnplus1(Hnplus1<0.02)=0.02;

% Calc Ns from soil production here then hand off to cosmo() after
% erosion(). 
Ns_soil_prod = 0; %((xr(host_min).*Ns*P.*exp(-alpha.*H))./(rhos.*X(:,:,host_min).*H))*dt + ((xr(host_min).*Nzb*P.*exp(-alpha.*H))./(rhos.*X(:,:,host_min).*H))*dt;

% boundary condition: no soil is produced, and therefore no bedrock
% lowering due to soil production, at boundaries.

Hnplus1(C==1) = H(C==1);
Hnplus1(C==2) = H(C==2);
Hnplus1(M2==1) = H(M2==1);
Hnplus1(M==1) = H(M==1);

%% Update bedrock elevation B
% Conservation of mass gives new bedrock surface elevation
B_minus1 = B;
B = B + (rhos/rhor)*(H - Hnplus1);

% Keep baselevel @ 0 m



BLV=(P/rhor)*exp(-alpha.*Hnplus1);
BLV(C==1)=E;
BLV(C==2)=0;
BLV(M==1)=E;

if sum(isnan(BLV(:)))>0
    fprintf('SoilProd.m: NaN(s) detected in BLV\n')
end

%% Update soil mineral abundances X
for i=1:nX
    
    Xtemp = X(:,:,i).*H./Hnplus1 + xr(i)*(Hnplus1-H)./Hnplus1;
    
    % Apply boundary conditions.
    Xtemp2 = X(:,:,i);
    Xtemp(C==1) = Xtemp2(C==1); % enforce boundary condition: no change
    Xtemp(M2==1) = Xtemp2(M2==1);
    Xtemp(M==1) = Xtemp(M==1);
    Xtemp(C==2) = Xtemp(C==2);
    % in X at boundaries
    X(:,:,i) = Xtemp;
    X(X<=0)=0;
end

if sum(isnan(X(:)))>0
    fprintf('SoilProd.m: NaN(s) detected in X\n')
end

%% Finish
% assign soil thickness and cosmo conc. at time t + dt to output argument
%H(M==0) = Hnplus1(M==0);
H(M2==0) = Hnplus1(M2==0);
H(C==1) = 1e-10;
H(M==1) = 1e-10;
%fprintf('SoilProd.m: Max X - %s\n', max(X(:)))
