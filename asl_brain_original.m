%% ASL_BRAIN - orignal version used for the PV correction paper

% load GM and WM maps
[gm_map,dims,scales] = ra('gm_asl');
wm_map = ra('wm_asl');
% load brain mask - gm+wm then fslmaths -bin
mask = ra('brain_mask');

gm_vec = vols2matrix(gm_map,mask);
wm_vec = vols2matrix(wm_map,mask);

snr = 15; %this is SNR for a curve that is for pure GM 60 ml/100g/min


%% perfusion parameters
cbf_gm = 60;
cbf_wm = 20;
tau_gm = 1.1;
tau_wm = 1.1;
at_gm = 0.7;
at_wm = 1;
T1_gm = 1.3;
T1_wm = 1.1;
T1b = 1.6;
lam_gm = 0.98;
lam_wm = 0.82;

gm_cbf_sd = 0;
wm_cbf_sd = 0;

gm_cbf = cbf_gm + gm_cbf_sd*randn(size(gm_vec));
wm_cbf = cbf_wm + wm_cbf_sd*randn(size(wm_vec));
gm_at = at_gm + 0*randn(size(gm_vec));
wm_at = at_wm + 0*randn(size(wm_vec));
gm_tau = tau_gm + 0*randn(size(gm_vec));
wm_tau= tau_wm + 0*randn(size(wm_vec));


%% curves
tis = 0.3:0.2:2.5;
src = zeros(length(gm_vec),length(tis));
% for time being loop
for i=1:length(gm_vec)
    kc_gm = eagle_box([gm_cbf(i),gm_at(i),gm_tau(i),0,0,0,0,T1_gm,T1b,lam_gm],tis,1,0);
    kc_wm = eagle_box([wm_cbf(i),wm_at(i),wm_tau(i),0,0,0,0,T1_gm,T1b,lam_wm],tis,1,0);
    src(i,:) = gm_vec(i)*kc_gm+wm_vec(i)*kc_wm;
end

%% add noise

%determine from the **flat** src data what the curve maximum is for pure GM
signalmag = max(kc_gm);
noise_sd = signalmag/snr
signal = src + noise_sd*randn(size(src));

%% make into brain
asl_brain = matrix2vols(signal,mask);
save_avw(asl_brain,'asl_brain','f',scales);

%%
gm_cbf_map = matrix2vols(gm_cbf,mask);
save_avw(gm_cbf_map,'asl_brain_flat_truegm','f',scales);
wm_cbf_map = matrix2vols(wm_cbf,mask);
save_avw(wm_cbf_map,'asl_brain_flat_truewm','f',scales);


%% add some spatial variation in cbf

x = (1:size(mask,1))/size(mask,1)*pi - pi/2;
y = (1:size(mask,2))/size(mask,2)*pi - pi/2;
z = (1:size(mask,3))/size(mask,3)*pi - pi/2;

gm_cbf_map = matrix2vols(gm_cbf,mask);
cbfsinx = repmat( 25*sin(x*4)' ,[1 64 22]);
cbfsiny = repmat( 10*sin(y*6) ,[64 1 22]);
cbfsinz = 0; %repmat( shiftdim( 15*sin(z*4),-1) ,[64 64 1]);
gm_cbf_map = gm_cbf_map + cbfsinx + cbfsiny + cbfsinz;
gm_cbf = vols2matrix(gm_cbf_map,mask);

%% curves
tis = 0.3:0.2:2.5;
src = zeros(length(gm_vec),length(tis));
% for time being loop
for i=1:length(gm_vec)
    kc_gm = eagle_box([gm_cbf(i),gm_at(i),gm_tau(i),0,0,0,0,T1_gm,T1b,lam_gm],tis,1,0);
    kc_wm = eagle_box([wm_cbf(i),wm_at(i),wm_tau(i),0,0,0,0,T1_gm,T1b,lam_wm],tis,1,0);
    src(i,:) = gm_vec(i)*kc_gm+wm_vec(i)*kc_wm;
end

%% add noise
%noise_sd = 1; ** defined from flat data!
signal = src + noise_sd*randn(size(src));

%% make into brain
asl_brain = matrix2vols(signal,mask);
save_avw(asl_brain,'asl_brain_sin','f',scales);

%%
gm_cbf_map = matrix2vols(gm_cbf,mask);
save_avw(gm_cbf_map,'asl_brain_sin_truegm','f',scales);
wm_cbf_map = matrix2vols(wm_cbf,mask);
save_avw(wm_cbf_map,'asl_brain_sin_truewm','f',scales);


%% add hyper/hypo-perfusion lesions
cbf_les = 20;
gm_cbf_map = matrix2vols(gm_cbf,mask);
[x,y,z] = meshgrid(1:size(mask,1),1:size(mask,2),1:size(mask,3));
cent = [21, 18, 11];
r = 10;
dr = sqrt( (x-cent(1)).^2 + (y-cent(2)).^2 + (z-cent(3)).^2);
[ind] = find(dr<=r);
gm_cbf_map(ind)=cbf_les;

% cbf_les = 100;
% cent = [42, 45, 11];
% r = 8;
% dr = sqrt( (x-cent(1)).^2 + (y-cent(2)).^2 + (z-cent(3)).^2);
% [ind] = find(dr<=r);
% gm_cbf_map(ind)=cbf_les;
% 
% need to get the altered CBF back into the vector we will be using next!
gm_cbf = vols2matrix(gm_cbf_map,mask);
%% curves
tis = 0.3:0.2:2.5;
src = zeros(length(gm_vec),length(tis));
% for time being loop
for i=1:length(gm_vec)
    kc_gm = eagle_box([gm_cbf(i),gm_at(i),gm_tau(i),0,0,0,0,T1_gm,T1b,lam_gm],tis,1,0);
    kc_wm = eagle_box([wm_cbf(i),wm_at(i),wm_tau(i),0,0,0,0,T1_gm,T1b,lam_wm],tis,1,0);
    src(i,:) = gm_vec(i)*kc_gm+wm_vec(i)*kc_wm;
end

%% add noise
%noise_sd = 1; **defined from flat GM
signal = src + noise_sd*randn(size(src));

%% make into brain
asl_brain = matrix2vols(signal,mask);
save_avw(asl_brain,'asl_brain_les','f',scales);

%%
gm_cmf_map = matrix2vols(gm_cbf,mask);
save_avw(gm_cbf_map,'asl_brain_les_truegm','f',scales);
wm_cbf_map = matrix2vols(wm_cbf,mask);
save_avw(wm_cbf_map,'asl_brain_les_truewm','f',scales);




