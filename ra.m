% This function reads in nifty file and determines the dimension
% Input:
% Full file name in string format
% Output:
% asl_4D_matrix: 4D matrix of ASL time series
% dims: Dimension of the 4D matrix
% scales: dummy variable, to be used in the future

function [asl_4D_matrix, dims, scales] = ra(full_file_name)

	scales = 0; % To be added

	file_handle = load_nii(full_file_name);
	asl_4D_matrix = file_handle.img;
	dims = size(input_4D_matrix);


end

