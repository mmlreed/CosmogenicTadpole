% NL_NPS nonlinear hillslope sediment transport nine point scheme with flux correction to 
% keep soil depth from going either negative or wildly high when slope values are close to
% critical. Flux is calculated as the sum of the eight neighbors of a cell (see
% Fagherazzi et al., 2002 WRR; Perron, 2011 JGR-ES). 

function [S, Q_k_tot, qns_tot, qx_tot] = NL_NPS_X_NS(S,C,M,M2,H,X,Ns,K,J,dt,dx,dy,S_c,K3,nX,rhos,host_min) 

%dt = 10;
%dx = 2;
%dy = 2;
del_x = dy;
del_y = dx;
diag_k = sqrt(dx^2 + dy^2);
diag_delk = (2*dx*dy)/diag_k;

% Critical slope - need to put this in with p struct
%S_c = 1.25; % Roering et al., 1999
% transport coefficient 
%K3 = 0.0032; % Roering et al., 1999

if sum(isnan(H(:)))>0
    fprintf('NL_NPS.m: NaN(s) detected in H before diffusion!\n')
end

if sum(isnan(S(:)))>0
    fprintf('NL_NPS.m: NaN(s) detected in S before diffusion!\n')
end

% save previous 
Snminus_1 = S;

% Let's call these right way from now on... (x,y) - (row,col)
% mirror y direction
% top side
S_y_L = [S(:,2) S];
% bottom side
S_y_R = [S S(:,end-1)];

% mirror x direction
% left side
S_x_T = [S(2,:); S]; % S(end-1,:)]; 
% right side
S_x_B = [S;  S(end-1,:)];


% mirrored for easier diagonals
S_diag = [S(:,2) S S(:,end-1)];
S_diag = [S_diag(2,:); S_diag; S_diag(end-1,:)]; 

S_xy_TL = S_diag(1:end-1, 1:end-1);
S_xy_BL = S_diag(2:end, 1:end-1);
S_xy_BR = S_diag(2:end, 2:end);
S_xy_TR = S_diag(1:end-1, 2:end);



% Need to compute spatial 1st derivatives for all eight directions with one-way
% differencing (z_ij - z_k)/dk going in counter-clockwise direction starting 
% from directly above z_ij. This choice is abritrary

S_k1 = (S_x_T(2:end, 1:end) - S_x_T(1:end-1, 1:end))/dx;

S_k2 = (S_xy_TL(2:end, 2:end) - S_xy_TL(1:end-1, 1:end-1))/diag_k; 

S_k3 = (S_y_L(1:end, 2:end) - S_y_L(1:end, 1:end-1))/dy;

S_k4 = (S_xy_BL(1:end-1,2:end) - S_xy_BL(2:end, 1:end-1))/diag_k;

S_k5 = (S_x_B(1:end-1, 1:end) - S_x_B(2:end, 1:end))/dx;

S_k6 = (S_xy_BR(1:end-1,1:end-1) - S_xy_BR(2:end, 2:end))/diag_k;

S_k7 = (S_y_R(1:end, 1:end-1) - S_y_R(1:end, 2:end))/dy;

S_k8 = (S_xy_TR(2:end, 1:end-1) - S_xy_TR(1:end-1,2:end))/diag_k;

% Slopes at 1.25 or above will produce infinite fluxes leading to NaNs 
% Need to may this generic

S_k1(S_k1>S_c)=S_c - 0.01;
S_k1(S_k1<-S_c)=-S_c + 0.01;

S_k2(S_k2>S_c)=S_c - 0.01;
S_k2(S_k2<-S_c)=-S_c + 0.01;

S_k3(S_k3>S_c)=S_c - 0.01;
S_k3(S_k3<-S_c)=-S_c + 0.01;

S_k4(S_k4>S_c)=S_c - 0.01;
S_k4(S_k4<-S_c)=-S_c + 0.01;

S_k5(S_k5>S_c)=S_c - 0.01;
S_k5(S_k5<-S_c)=-S_c + 0.01;

S_k6(S_k6>S_c)=S_c - 0.01;
S_k6(S_k6<-S_c)=-S_c + 0.01;

S_k7(S_k7>S_c)=S_c - 0.01;
S_k7(S_k7<-S_c)=-S_c + 0.01;

S_k8(S_k8>S_c)=S_c - 0.01;
S_k8(S_k8<-S_c)=-S_c + 0.01;

% calculate vector flux qk for each direction using nonlinear diffusion

q_k1 = (-K3*S_k1)./(1-(S_k1/S_c).^2);

q_k2 = (-K3*S_k2)./(1-(S_k2/S_c).^2);

q_k3 = (-K3*S_k3)./(1-(S_k3/S_c).^2);

q_k4 = (-K3*S_k4)./(1-(S_k4/S_c).^2);

q_k5 = (-K3*S_k5)./(1-(S_k5/S_c).^2);

q_k6 = (-K3*S_k6)./(1-(S_k6/S_c).^2);

q_k7 = (-K3*S_k7)./(1-(S_k7/S_c).^2);

q_k8 = (-K3*S_k8)./(1-(S_k8/S_c).^2);

% No flux in streams

q_k1(C==1)=0;
q_k2(C==1)=0;
q_k3(C==1)=0;
q_k4(C==1)=0;
q_k5(C==1)=0;
q_k6(C==1)=0;
q_k7(C==1)=0;
q_k8(C==1)=0;

% No flux in lake 

q_k1(C==2)=0;
q_k2(C==2)=0;
q_k3(C==2)=0;
q_k4(C==2)=0;
q_k5(C==2)=0;
q_k6(C==2)=0;
q_k7(C==2)=0;
q_k8(C==2)=0;

% No flux at edges -- may be redundant with below

q_k1(M2==1)=0;
q_k2(M2==1)=0;
q_k3(M2==1)=0;
q_k4(M2==1)=0;
q_k5(M2==1)=0;
q_k6(M2==1)=0;
q_k7(M2==1)=0;
q_k8(M2==1)=0;

% No flux outside of watershed

q_k1(M==1)=0;
q_k2(M==1)=0;
q_k3(M==1)=0;
q_k4(M==1)=0;
q_k5(M==1)=0;
q_k6(M==1)=0;
q_k7(M==1)=0;
q_k8(M==1)=0;




% turn these into Qk's; Qk = del_k*qk

Q_k1 = q_k1*dy;

Q_k2 = q_k2*diag_delk;

Q_k3 = q_k3*dx;

Q_k4 = q_k4*diag_delk;

Q_k5 = q_k5*dy;

Q_k6 = q_k6*diag_delk;

Q_k7 = q_k7*dx;

Q_k8 = q_k8*diag_delk;

dz_dt = (1/(dx*dy))*(Q_k1 + Q_k2 + Q_k3 + Q_k4 + Q_k5 + Q_k6 + Q_k7 + Q_k8)*dt;

dz_dt(C==1)=0;
dz_dt(C==2)=0;
dz_dt(M==1)=0;
dz_dt(M2==1)=0;

dz_dt_orig = dz_dt;


if sum(isnan(dz_dt(:)))>0
    fprintf('NL_NPS.m: NaN(s) detected in dz_dt\n')
end

% Correct for high outgoing fluxes -- This does "lock in" very low soil
% depths in areas with slopes near the critical slope (a feature not a
% bug?).. The mininum soil thickness of 0.01 can be changed... this does get checked elsewhere to
% prevent some divide by zero operations

qs_ratio=ones(size(S));

for n=1:K
    for m=1:J
        if dz_dt(n,m)<-1*H(n,m)
            dz_dt(n,m)=-1*H(n,m);
            qs_ratio(n,m)=dz_dt_orig(n,m)/dz_dt(n,m);
        end
    end
end

qs_list=find(qs_ratio>0);
 
% Need to optimize here
for i=1:length(qs_list)
    idx = qs_list(i);
    if q_k1(idx)<0
        q_k1(idx)=q_k1(idx)/qs_ratio(idx);
    end
end

for i=1:length(qs_list)
    idx = qs_list(i);
    if q_k2(idx)<0
        q_k2(idx)=q_k2(idx)/qs_ratio(idx);
    end
end

for i=1:length(qs_list)
    idx = qs_list(i);
    if q_k3(idx)<0
        q_k3(idx)=q_k3(idx)/qs_ratio(idx);
    end
end

for i=1:length(qs_list)
    idx = qs_list(i);
    if q_k4(idx)<0
        q_k4(idx)=q_k4(idx)/qs_ratio(idx);
    end
end

for i=1:length(qs_list)
    idx = qs_list(i);
    if q_k5(idx)<0
        q_k5(idx)=q_k5(idx)/qs_ratio(idx);
    end
end

for i=1:length(qs_list)
    idx = qs_list(i);
    if q_k6(idx)<0
        q_k6(idx)=q_k6(idx)/qs_ratio(idx);
    end
end

for i=1:length(qs_list)
    idx = qs_list(i);
    if q_k7(idx)<0
        q_k7(idx)=q_k7(idx)/qs_ratio(idx);
    end
end

for i=1:length(qs_list)
    idx = qs_list(i);
    if q_k8(idx)<0
        q_k8(idx)=q_k8(idx)/qs_ratio(idx);
    end
end

% Now deposition fluxes must be recalculated
% Loop through everything again to correct for high incoming fluxes; We
% don't care about flux correction at locations where M2==1, so these
% points are ignored

for n=2:K-1
    for m=2:J-1
        if q_k1(n,m)~=-q_k5(n-1,m) & q_k1(n,m)>0
           q_k1(n,m)=-q_k5(n-1,m);
        end
    end
end

for n=2:K-1
    for m=2:J-1
        if q_k2(n,m)~=-q_k6(n-1,m-1) & q_k2(n,m)>0
           q_k2(n,m)=-q_k6(n-1,m-1);
        end
    end
end

for n=2:K-1
    for m=2:J-1
        if q_k3(n,m)~=-q_k7(n,m-1) & q_k3(n,m)>0
           q_k3(n,m)=-q_k7(n,m-1);
        end
    end
end

for n=2:K-1
    for m=2:J-1
        if q_k4(n,m)~=-q_k8(n+1,m-1) & q_k4(n,m)>0
           q_k4(n,m)=-q_k8(n+1,m-1);
        end
    end
end

for n=2:K-1
    for m=2:J-1
        if q_k5(n,m)~=-q_k1(n+1,m) & q_k5(n,m)>0
           q_k5(n,m)=-q_k1(n+1,m);
        end
    end
end

for n=2:K-1
    for m=2:J-1
        if q_k6(n,m)~=-q_k2(n+1,m+1) & q_k6(n,m)>0
           q_k6(n,m)=-q_k2(n+1,m+1);
        end
    end
end

for n=2:K-1
    for m=2:J-1
        if q_k7(n,m)~=-q_k3(n,m+1) & q_k7(n,m)>0
           q_k7(n,m)=-q_k3(n,m+1);
        end
    end
end

for n=2:K-1
    for m=2:J-1
        if q_k8(n,m)~=-q_k4(n-1,m+1) & q_k8(n,m)>0
           q_k8(n,m)=-q_k4(n-1,m+1);
        end
    end
end


Q_k1 = q_k1*dy;

Q_k2 = q_k2*diag_delk;

Q_k3 = q_k3*dx;

Q_k4 = q_k4*diag_delk;

Q_k5 = q_k5*dy;

Q_k6 = q_k6*diag_delk;

Q_k7 = q_k7*dx;

Q_k8 = q_k8*diag_delk;

Q_k_tot = (Q_k1 + Q_k2 + Q_k3 + Q_k4 + Q_k5 + Q_k6 + Q_k7 + Q_k8);

dz_dt = (1/(dx*dy))*(Q_k1 + Q_k2 + Q_k3 + Q_k4 + Q_k5 + Q_k6 + Q_k7 + Q_k8)*dt;


% No change zones


dz_dt(C==1)=0;

dz_dt(C==2)=0;

dz_dt(M2==1)=0;

dz_dt(M==1)=0;

S = Snminus_1 + dz_dt;


if sum(isnan(S(:)))>0
    fprintf('NL_NPS.m: NaN(s) detected in S\n')
end

% For mineral and cosmo transport/conc. conservation

qs_tot = (q_k1 + q_k2 + q_k3 + q_k4 + q_k5 + q_k6 + q_k7 + q_k8);

% Calc changes in X due to transport in terms of mass (+ -> gain; - ->
% loss)

qx_k1=zeros(size(X));
qx_k2=zeros(size(X));
qx_k3=zeros(size(X));
qx_k4=zeros(size(X));
qx_k5=zeros(size(X));
qx_k6=zeros(size(X));
qx_k7=zeros(size(X));
qx_k8=zeros(size(X));

for i=1:nX
qx_k1(:,:,i) = q_k1.*dx.*rhos.*X(:,:,i);
qx_k2(:,:,i) = q_k2.*diag_k.*rhos.*X(:,:,i);
qx_k3(:,:,i) = q_k3.*dy.*rhos.*X(:,:,i);
qx_k4(:,:,i) = q_k4.*diag_k.*rhos.*X(:,:,i);
qx_k5(:,:,i) = q_k5.*dx.*rhos.*X(:,:,i);
qx_k6(:,:,i) = q_k6.*diag_k.*rhos.*X(:,:,i);
qx_k7(:,:,i) = q_k7.*dy.*rhos.*X(:,:,i);
qx_k8(:,:,i) = q_k8.*diag_k.*rhos.*X(:,:,i);
end

% positive fluxes require usage of data from point of origin instead of local; allocate these at there are threeeee lewps

qx_k1_p=zeros(size(X));
qx_k2_p=zeros(size(X));
qx_k3_p=zeros(size(X));
qx_k4_p=zeros(size(X));
qx_k5_p=zeros(size(X));
qx_k6_p=zeros(size(X));
qx_k7_p=zeros(size(X));
qx_k8_p=zeros(size(X));

% Very costly loops here..

for k=2:K-1
    for j=2:J-1
        for i=1:nX
            qx_k1_p(k,j,i) = q_k1(k,j).*dx.*rhos.*X(k-1,j,i);
            qx_k2_p(k,j,i) = q_k2(k,j).*diag_k.*rhos.*X(k-1,j-1,i);
            qx_k3_p(k,j,i) = q_k3(k,j).*dy.*rhos.*X(k,j-1,i);
            qx_k4_p(k,j,i) = q_k4(k,j).*diag_k.*rhos.*X(k+1,j-1,i);
            qx_k5_p(k,j,i) = q_k5(k,j).*dx.*rhos.*X(k+1,j,i);
            qx_k6_p(k,j,i) = q_k6(k,j).*diag_k.*rhos.*X(k+1,j+1,i);
            qx_k7_p(k,j,i) = q_k7(k,j).*dy.*rhos.*X(k,j+1,i);
            qx_k8_p(k,j,i) = q_k8(k,j).*diag_k.*rhos.*X(k-1,j+1,i);
        end
    end
end

qx_k1(q_k1>0)=qx_k1_p(q_k1>0);
qx_k2(q_k2>0)=qx_k2_p(q_k2>0);
qx_k3(q_k3>0)=qx_k3_p(q_k3>0);
qx_k4(q_k4>0)=qx_k4_p(q_k4>0);
qx_k5(q_k5>0)=qx_k5_p(q_k5>0);
qx_k6(q_k6>0)=qx_k6_p(q_k6>0);
qx_k7(q_k7>0)=qx_k7_p(q_k7>0);
qx_k8(q_k8>0)=qx_k8_p(q_k8>0);


% Total mass flux in minerals
qx_tot = (qx_k1 + qx_k2 + qx_k3 + qx_k4 + qx_k5 + qx_k6 + qx_k7 + qx_k8);

% No flux zones
qx_tot(C==1)=0;
qx_tot(C==2)=0;
qx_tot(M2==1)=0;
qx_tot(M==1)=0;


% Flux of 10Be atoms in all 8 direction; Quartz is fluxed using the same
% volumetric flux, so these can be modeled separately and still be
% conserved

qns_k1 = q_k1.*dx.*rhos.*X(:,:,host_min).*Ns;
qns_k2 = q_k2.*diag_k.*X(:,:,host_min).*rhos.*Ns;
qns_k3 = q_k3.*dy.*X(:,:,host_min).*rhos.*Ns;
qns_k4 = q_k4.*diag_k.*rhos.*X(:,:,host_min).*Ns;
qns_k5 = q_k5.*dx.*rhos.*X(:,:,host_min).*Ns;
qns_k6 = q_k6.*diag_k.*rhos.*X(:,:,host_min).*Ns;
qns_k7 = q_k7.*dy.*rhos.*X(:,:,host_min).*Ns;
qns_k8 = q_k8.*diag_k.*X(:,:,host_min).*rhos.*Ns;

% allocate
qns_k1_p=zeros(size(C));
qns_k2_p=zeros(size(C));
qns_k3_p=zeros(size(C));
qns_k4_p=zeros(size(C));
qns_k5_p=zeros(size(C));
qns_k6_p=zeros(size(C));
qns_k7_p=zeros(size(C));
qns_k8_p=zeros(size(C));

% Correct for positive fluxes into grid cell

for k=2:K-1
    for j=2:J-1
            qns_k1_p(k,j) = q_k1(k,j).*dx.*rhos.*X(k-1,j,host_min).*Ns(k-1,j);
            qns_k2_p(k,j) = q_k2(k,j).*diag_k.*rhos.*X(k-1,j-1,host_min).*Ns(k-1,j-1);
            qns_k3_p(k,j) = q_k3(k,j).*dy.*rhos.*X(k,j-1,host_min).*Ns(k,j-1);
            qns_k4_p(k,j) = q_k4(k,j).*diag_k.*rhos.*X(k+1,j-1,host_min).*Ns(k+1,j-1);
            qns_k5_p(k,j) = q_k5(k,j).*dx.*rhos.*X(k+1,j,host_min).*Ns(k+1,j);
            qns_k6_p(k,j) = q_k6(k,j).*diag_k.*rhos.*X(k+1,j+1,host_min).*Ns(k+1,j+1);
            qns_k7_p(k,j) = q_k7(k,j).*dy.*rhos.*X(k,j+1,host_min).*Ns(k,j+1);
            qns_k8_p(k,j) = q_k8(k,j).*diag_k.*rhos.*X(k-1,j+1,host_min).*Ns(k-1,j+1);
    end
end


qns_k1(q_k1>0)=qns_k1_p(q_k1>0);
qns_k2(q_k2>0)=qns_k2_p(q_k2>0);
qns_k3(q_k3>0)=qns_k3_p(q_k3>0);
qns_k4(q_k4>0)=qns_k4_p(q_k4>0);
qns_k5(q_k5>0)=qns_k5_p(q_k5>0);
qns_k6(q_k6>0)=qns_k6_p(q_k6>0);
qns_k7(q_k7>0)=qns_k7_p(q_k7>0);
qns_k8(q_k8>0)=qns_k8_p(q_k8>0);

qns_tot = (qns_k1 + qns_k2 + qns_k3 + qns_k4 + qns_k5 + qns_k6 + qns_k7 + qns_k8);

% No flux zones
qns_tot(C==1)=0;
qns_tot(C==2)=0;
qns_tot(M2==1)=0;
qns_tot(M==1)=0;


end
































