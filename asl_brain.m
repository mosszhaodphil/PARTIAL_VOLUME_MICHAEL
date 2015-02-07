%% ASL_BRAIN

srcdir='brain_src/'

% load PVE maps (in ASL data space)
[gm_map,dims,scales] = ra([srcdir 'gm_asl']);
wm_map = ra([srcdir 'wm_asl']);
csf_map = ra([srcdir 'csf_asl']);

% load brain mask
mask = ra([srcdir 'mask']);

gm_vec = vols2matrix(gm_map,mask);
wm_vec = vols2matrix(wm_map,mask);
csf_vec = vols2matrix(csf_map,mask);


%% Parameters

% perfusion
cbf_gm = 60 /6000; % s^-1
cbf_wm = 20 /6000; % s^-1
tau_gm = 1.1; % s
tau_wm = 1.1; % s
at_gm = 0.7; % s
at_wm = 1; % s

% constants
T1_gm = 1.3;
T1_wm = 1.1;
T1_csf = 4.3;
T1b = 1.6;
lam_gm = 0.98;
lam_wm = 0.82;
lam_csf = 0.9;

% static
M0b = 10000;

%% [OPTION] random voxelwise variations in perfusion parameters
gm_cbf_sd = 0;
wm_cbf_sd = 0;

gm_cbf = cbf_gm + gm_cbf_sd*randn(size(gm_vec));
wm_cbf = cbf_wm + wm_cbf_sd*randn(size(wm_vec));
gm_at = at_gm + 0*randn(size(gm_vec));
wm_at = at_wm + 0*randn(size(wm_vec));
gm_tau = tau_gm + 0*randn(size(gm_vec));
wm_tau= tau_wm + 0*randn(size(wm_vec));

%% [OPTION] add some spatial variation in cbf

x = (1:size(mask,1))/size(mask,1)*pi - pi/2;
y = (1:size(mask,2))/size(mask,2)*pi - pi/2;
z = (1:size(mask,3))/size(mask,3)*pi - pi/2;

gm_cbf_map = matrix2vols(gm_cbf,mask);
cbfsinx = repmat( 25*sin(x*4)' ,[1 64 22]);
cbfsiny = repmat( 10*sin(y*6) ,[64 1 22]);
cbfsinz = 0; %repmat( shiftdim( 15*sin(z*4),-1) ,[64 64 1]);
gm_cbf_map = gm_cbf_map + cbfsinx + cbfsiny + cbfsinz;
gm_cbf = vols2matrix(gm_cbf_map,mask);

%% [OPTION] add hyper/hypo-perfusion lesions
cbf_les = 20;
gm_cbf_map = matrix2vols(gm_cbf,mask);
[x,y,z] = meshgrid(1:size(mask,1),1:size(mask,2),1:size(mask,3));
cent = [21, 18, 11];
r = 10;
dr = sqrt( (x-cent(1)).^2 + (y-cent(2)).^2 + (z-cent(3)).^2);
[ind] = find(dr<=r);
gm_cbf_map(ind)=cbf_les;

cbf_les = 100;
cent = [42, 45, 11];
r = 8;
dr = sqrt( (x-cent(1)).^2 + (y-cent(2)).^2 + (z-cent(3)).^2);
[ind] = find(dr<=r);
gm_cbf_map(ind)=cbf_les;

% need to get the altered CBF back into the vector we will be using next!
gm_cbf = vols2matrix(gm_cbf_map,mask);


%% kinetic curves
tis = 0.3:0.2:2.5;
src = zeros(length(gm_vec),length(tis));
% for time being loop
for i=1:length(gm_vec)
    kc_gm(i,:) = eagle_box([gm_cbf(i),gm_at(i),gm_tau(i),0,0,0,0,T1_gm,T1b,lam_gm],tis,1,0);
    kc_wm(i,:) = eagle_box([wm_cbf(i),wm_at(i),wm_tau(i),0,0,0,0,T1_gm,T1b,lam_wm],tis,1,0);
    %src(i,:) = gm_vec(i)*kc_gm+wm_vec(i)*kc_wm;
end

%% static tissue curves

% M0 for different tissue
M0gm = M0b*lam_gm;
M0wm = M0b*lam_wm;
M0csf = M0b*lam_csf;
A=0.9; %saturation efficiency

% control images
% saturation recovery
for i=1:length(gm_vec)
    ctrl(i,:) = M0gm*gm_vec(i)*(1 - A*exp(-tis/T1_gm)) + M0wm*wm_vec(i)*(1 - A*exp(-tis/T1_wm)) + M0csf*csf_vec(i)*(1 - A*exp(-tis/T1_csf));
end

% tag images
for i=1:length(gm_vec)
    tag(i,:) = ctrl(i,:) -  M0gm*gm_vec(i)*kc_gm(i,:) - M0wm*wm_vec(i)*kc_wm(i,:);
end

% combine tag and control
src(:,1:2:2*length(tis)) = tag;
src(:,2:2:2*length(tis)) = ctrl;

%% create dataset containing repeats
nrpts=5;
src = repmat(src,[1,nrpts]);

%% add noise - relative to M0b

snr = 1000;

noise_sd = M0b/snr;
signal = src + noise_sd*randn(size(src));

%% make into brain
asl_brain = matrix2vols(signal,mask);

%% save
save_avw(asl_brain,'asl_brain','f',scales);

% true CBF values
gm_cbf_map = matrix2vols(gm_cbf,mask);
save_avw(gm_cbf_map,'asl_brain_truegm','f',scales);
wm_cbf_map = matrix2vols(wm_cbf,mask);
save_avw(wm_cbf_map,'asl_brain_truewm','f',scales);


