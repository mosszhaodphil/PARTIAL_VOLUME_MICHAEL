
% Output file names
file_name_extracted_brain = 'str_brain';
file_name_registration = 'str_to_t1';
file_name_seg_gm = 'str_to_t1_pve_1';
file_name_seg_wm = 'str_to_t1_pve_2';
file_name_gm_raw = 'pvgm_raw';
file_name_wm_raw = 'pvwm_raw';
file_name_gm = 'pvgm';
file_name_wm = 'pvwm';
file_name_gm_mask = 'gm_mask';
file_name_wm_mask = 'wm_mask';


% Extract brain
call_fsl(strcat('bet ', str_file, ' ', file_name_extracted_brain))

% Rigis registration
call_fsl(strcat('flirt -dof 6 ', '-in ', file_name_extracted_brain, ' -ref ', t1_file, ' -out ', file_name_registration))

% Segmentation
call_fsl(strcat('fast -p ', file_name_registration))

% Resample
call_fsl(strcat('applywarp -i ', file_name_seg_gm, ' -o ', file_name_gm_raw, ' -s --interp=trilinear'))
call_fsl(strcat('applywarp -i ', file_name_seg_wm, ' -o ', file_name_wm_raw, ' -s --interp=trilinear'))

% Thresholding at 0.1
call_fsl(strcat('fslmaths ', file_name_gm_raw, ' -thr 0.1 -min 1 ', file_name_gm))
call_fsl(strcat('fslmaths ', file_name_wm_raw, ' -thr 0.1 -min 1 ', file_name_wm))

% Get binary tissue map
call_fsl(strcat('fslmaths ', file_name_gm, ' -bin ', file_name_gm_mask))
call_fsl(strcat('fslmaths ', file_name_wm, ' -bin ', file_name_wm_mask))

% Now you have generated 


