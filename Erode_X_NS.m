function [H, B, X, Ns] = Erode_X_NS(H,B,X,Ns,Nzb,C,nX,K,J,dt,dx,dy,S_c,K3,xr,M,M2,rhos,host_min)
%
% Update soil thickness H, bedrock elevation B, and mineral abundances X
% according to the perturbations to these quantities by soil erosion.
% This is part of a splitting method, whereby this function computes a
% portion of the temporal changes in H, B, and X, and other functions
% compute the other portions of those changes.


%% Update soil thickness H

% NONLINEAR DIFFUSION
S = B + H; % calculate surface elevation
%S = ADI(S,K,J,Aplus,Aminus,Bplus,Bminus,M2); % step surface elevation forward 
% in time due to hillslope soil transport
[S, Q_k_tot, qns_tot, qx_tot] = NL_NPS_X_NS(S,C,M,M2,H,X,Ns,K,J,dt,dx,dy,S_c,K3,nX,rhos,host_min);

Htemp=H;

if sum(isnan(S(:)))>0
    fprintf('Erode.m: NaN(s) detected in S\n')
end

Hnplus1 = S - B; % calculate new soil thickness
%Hnplus1(M==1)=Htemp(M==1);

if sum(isnan(H(:)))>0
    fprintf('Erode.m: NaN(s) detected in H\n')
end

Hnplus1(Hnplus1<0.02)=0.02;
Hnplus1(M2==1)=1e-10;
Hnplus1(M==1)=1e-10;
Hnplus1(C==1)=1e-10;
Hnplus1(C==2)=1e-10;


%% Update mineral abundances X

% Use vector flux mineral flux to calculate change in concentration

%dX_dt=qx_tot./(H*dx*dy*rhos);

X_old=X;

%Xnplus1 = X + dX_dt*dt;

Xnplus1 =  ((H*dx*dy*rhos.*X_old + qx_tot*dt)./(H*dx*dy*rhos + Q_k_tot*rhos*dt));

    
% Apply boundary conditions.
for i = 1:nX
    Xtemp = Xnplus1(:,:,i);
    Xtemp2 = X(:,:,i);
    Xtemp(C==1) = Xtemp2(C==1);
    Xtemp(C==2) = Xtemp2(C==2);
    % enforce boundary condition: no change
    Xtemp(M2==1) = Xtemp(M2==1);
    Xtemp(M==1) = Xtemp(M==1);
    % in X at boundaries.
    Xnplus1(:,:,i) = Xtemp;
end

%% Finish

if sum(isnan(X(:)))>0
    fprintf('Erode.m: NaN(s) detected in X\n')
end


H_before=H;

% assign soil thickness and mineral concentrations at time t + dt to output
% arguments
% H before diffusion as the old H is passed to ForwardEuler/RK4 we need to
% eval Ns using the pre-diffusion H to lessen error between instaneous and
% inferred rates
H = Hnplus1;
Xnplus1(Xnplus1<=0)=0;
X = Xnplus1;
X_ratio=ones(size(B));

% Percentage concentrations need to sum to 1; possibly vectorized?

for i=1:K
    for j=1:J
       if sum(X(i,j,:))>1
            X_ratio(i,j)=sum(X(i,j,:));
            X(i,j,:)=X(i,j,:)/X_ratio(i,j);
        end
    end
end

for i=1:K
    for j=1:J
        
        if sum(X(i,j,:))<1
            X_ratio(i,j)=sum(X(i,j,:));
            X(i,j,:)=X(i,j,:)/X_ratio(i,j);
        end
    end
end

% for i=1:K
%     for j=1:J
%        if sum(X_alt(i,j,:))>1
%             X_alt_ratio(i,j)=sum(X_alt(i,j,:));
%             X_alt(i,j,:)=X_alt(i,j,:)/X_alt_ratio(i,j);
%         end
%     end
% end
% 
% for i=1:K
%     for j=1:J
%         if sum(X_alt(i,j,:))<1
%             X_alt_ratio(i,j)=sum(X_alt(i,j,:));
%             X_alt(i,j,:)=X_alt(i,j,:)/X_alt_ratio(i,j);
%         end
%     end
% end


Ns = ((Ns.*H_before*dx*dy*rhos.*X_old(:,:,host_min) + qns_tot*dt)./(H_before*dx*dy*rhos.*X_old(:,:,host_min) + Q_k_tot*rhos.*X(:,:,host_min)*dt));

% Set Ns to bedrock if H dips into minimum (think a thin mantle of gravel)
Ns(H==0.02)=Nzb(H==0.02);
Ns(C==1)=0;
Ns(C==2)=0;
Ns(M2==1)=0;
Ns(Ns<0)=0;

% This can be vectorized ; looks redudant with above

for i=1:K
    for j=1:J
        if H(i,j)==0.02
            X(i,j,:)=xr;
            Ns(i,j)=Nzb(i,j);
        end
    end
end

%fprintf('erode - soil 10Be concentration at (200,118): %s\n', Ns(200,118))

% fprintf('Erode.m: Max X - %s\n', max(X(:)))
% fprintf('Erode.m: Ns(50,50) - %s\n', Ns(50,50))
%fprintf('Erode.m: Ns_alt2(50,50) - %s\n', Ns_alt2(50,50))

