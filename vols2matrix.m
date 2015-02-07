% This function converts a 4D matrix into 2D matrix with mask
% The 4D matrix is read from nifty file
% The input 3D define the dimension of the output matrix.
% The input 4D matrix projects values to the 3D mask.
% The output is a 2D matrix of time series on each row.
% Input:
% input_4D_matrix: input 4D matrix read from nifty file
% input_3D_mask: user define mask
% output_matrix: output 2D matrix with ASl time series signal on each row

function output_matrix = vols2matrix(input_4D_matrix, input_3D_mask)

	% x y z dimension of input_4D_matrix must be the same with the x y z dimension of input_3D_matrix
	
	% process to check dimension of input matrices
	% To be implemented

	[x, y, z, t] = size(input_4D_matrix); % get dimension of 3D mask to be used in output matrix dimension

	output_matrix = zeros(x * y * z, t);

	counter = 1;

	for i = 1 : x
		for j = 1 : y
			for k = 1 : z
				output_matrix(counter, :) = input_4D_matrix(i, j, k, :);
				counter = counter + 1;
			end
		end
	end

end

