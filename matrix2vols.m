% This function converts a 2D matrix into 4D matrix based on the 3D mask file
% The 2D matrix contains the time series of ASL signal
% The output 4D matrix will be used to construct nifty file
% Input:
% input_matrix: input 2D matrix
% input_3D_mask: input 3D mask
% Output:
% output_4D_matrix: 4D matrix to be used to construct nifty file

function output_4D_matrix = matrix2vols(input_matrix, input_3D_mask)

	[x, y, z] = size(input_3D_mask); % get dimension of 3D mask matrix 
	[m, n] = size(input_matrix); % get dimension of 2D input matrix

	output_4D_matrix = zeros(x, y, z, n); % create an empty whose dimension is x y z n. n is the number of sampling points

	counter = 1;
	for i = 1 : x
		for j = 1 : y
			for k = 1 : z
				output_4D_matrix(i, j, k, :) = input_matrix(counter, :); % Assign the current row to output matrix's voxel
				counter = counter + 1;
			end
		end
	end

end

