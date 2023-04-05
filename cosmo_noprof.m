% Cosmo operations where getting pretty lengthy so I am making a separate
% m-file....This is not part of the previous operator-splitting scheme, but
% its own forward scheme with various sub-schemes(?). Still needs chemical
% erosion and decay. 


function [Ns, Nzb, Nzb_spal, Nzb_neg, Nzb_fast, Dinf_m, Dinf_m_lal, Dinf_m_zb] = cosmo_noprof(H, Ns, Nzb, Nzb_spal, Nzb_neg, Nzb_fast, C, rhor, rhos, cosmo_prod_spal, cosmo_prod_neg_1,...
    cosmo_prod_fast, lambda_10Be, L1, L2, L4, dt, t, K, J, xr, X, P, alpha, M2, M, host_min, zr_min)


Ns_spal = ((cosmo_prod_spal.*L1.*(1-exp(-rhos.*H./L1))));
Ns_neg_1 = ((cosmo_prod_neg_1.*L2.*(1-exp(-rhos.*H./L2)))); 
% Ns_neg_2 = ((cosmo_prod_neg_2*L3.*(1-exp(-rhos.*H/L3))));
Ns_fast = ((cosmo_prod_fast.*L4.*(1-exp(-rhos.*H./L4))));
Ns_prod_tot = (Ns_spal + Ns_neg_1 + Ns_fast)./(rhos.*H); % If using two expoentials for neg muons add here

Ns_1 = Ns + Ns_prod_tot*dt - ((xr(host_min).*Ns*P.*exp(-alpha.*H))./(rhos.*X(:,:,host_min).*H))*dt + ((xr(host_min).*Nzb*P.*exp(-alpha.*H))./(rhos.*X(:,:,host_min).*H))*dt - lambda_10Be*Ns*dt;
% %Ns_1 = Ns + Ns_prod_tot*dt - Ns_soil_prod - lambda_10Be*Ns*dt;
% 
Ns(M==0)=Ns_1(M==0); % Everything outside of catchment remains the same

Nzb_prod_spal = cosmo_prod_spal.*exp(-rhos.*H./L1)*dt;
Nzb_prod_neg =  cosmo_prod_neg_1.*exp(-rhos.*H./L2)*dt;
Nzb_prod_fast = cosmo_prod_fast.*exp(-rhos.*H./L4)*dt;


Nzb_spal = Nzb_spal + Nzb_prod_spal - P.*exp(-alpha.*H).*(Nzb_spal./L1)*dt - Nzb_spal*lambda_10Be*dt;
Nzb_neg = Nzb_neg + Nzb_prod_neg - P.*exp(-alpha.*H).*(Nzb_neg./L2)*dt - Nzb_neg*lambda_10Be*dt;
Nzb_fast = Nzb_fast + Nzb_prod_fast - P.*exp(-alpha.*H).*(Nzb_fast./L4)*dt - Nzb_fast*lambda_10Be*dt;



Nzb = Nzb_spal + Nzb_neg + Nzb_fast;

 
% Zero out soil in channels, boundaries, and lake cells 

Ns(M2==1)=0;
Ns(C==1) = 0;
Ns(C==2) = 0;
Ns(M==1)=0;

 
% % Zero out bedrock in channels, boundaries, and lake cells

Nzb(M2==1)=0;
Nzb(C==1)=0;
Nzb(C==2)=0;
Nzb(M==1)=0;

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