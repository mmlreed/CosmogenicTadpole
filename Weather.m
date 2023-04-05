function [H, X] = Weather(H,X,C,nX,kx,Ax,sx,wx,rhos,xr,dt,M2,M)
%
% Update soil thickness H and mineral abundances X according to the 
% perturbations to these quantities by mineral weathering. This is part of
% a splitting method, whereby this function computes a portion of the 
% temporal changes in H and X, and other functions compute the other 
% portions of those changes.


%% Update mineral abundances X

Xnplus1 = zeros(size(X));

if sum(isnan(X(:)))>0
    fprintf('Weather.m: NaN(s) detected in X at outset!\n')
end

% calculate the summation term
BigSigma = 0; 
for i=1:nX
    BigSigma = BigSigma + (kx(i)*Ax(i)*X(:,:,i) - sx(i)*wx(i)/rhos);
end

if sum(isnan(BigSigma(:)))>0
    fprintf('Weather.m: NaN(s) detected in BigSigma! [1st instance]\n')
end

% calculate Xn+1 for each mineral species
    
for i=1:nX
    Xtemp2 = X(:,:,i);
    Xtemp = X(:,:,i);
    DeltaX = dt * (-kx(i)*Ax(i)*X(:,:,i) + sx(i)*wx(i)/rhos + ...
        X(:,:,i).*BigSigma);
    Xtemp = Xtemp + DeltaX;

    % Apply boundary conditions.
    Xtemp(C==1) = Xtemp2(C==1);% enforce boundary condition: no change
    %Xtemp(M==1) = Xtemp2(M==1);
    Xtemp(M2==1) = Xtemp2(M2==1);
    Xtemp(C==2) = Xtemp2(C==2);
    % in X at boundaries

    % Finalize X.
    Xnplus1(:,:,i) = Xtemp;
end

if sum(isnan(Xnplus1(:)))>0
    fprintf('Weather.m: NaN(s) detected in Xnplus1!\n')
end

%% Update soil thickness H

% calculate the summation term at time n+1/2
BigSigma = 0; 
for i=1:nX
    BigSigma = BigSigma + (kx(i)*Ax(i)*0.5*(X(:,:,i)+Xnplus1(:,:,i)) - ...
        sx(i)*wx(i)/rhos);
end

if sum(isnan(BigSigma(:)))>0
    fprintf('Weather.m: NaN(s) detected in BigSigma! [2nd instance]\n')
end


Hnplus1 = H.*exp(-BigSigma*dt);

if sum(isnan(Hnplus1(:)))>0
    fprintf('Weather.m: NaN(s) detected in Hnplus1!\n')
end

Hnplus1(Hnplus1<0.02)=0.02;
% Boundary condition: no change in H due to weathering at boundaries.
Hnplus1(C==1) = H(C==1);
Hnplus1(C==2) = H(C==2); % ALL GOOD
Hnplus1(M==1) = H(M==1);
Hnplus1(M2==1) = H(M2==1);


%% Finish
% assign soil thickness and mineral concentrations at time t + dt to output
% arguments
H = Hnplus1;
Xnplus1(Xnplus1<=0)=0;
X = Xnplus1;
if sum(isnan(H(:)))>0
    fprintf('Weather.m: NaN(s) detected in H at end of Weather.m\n')
end
%fprintf('Weather.m: Max X - %s\n', max(X(:)))
