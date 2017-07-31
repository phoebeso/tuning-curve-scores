function correlationMatrix = calculate_correlation_matrix(matrix1,matrix2,n)
% Input:mat1, mat2 assumed to be the same size.
% n: matrix holding the number of cells joining multiplication in both matrices per time unit.
% output: correlation matrix
% If mat2 is not given it is assumed to be like mat1
% if n is not given it is calculated as the maximum possible

if (nargin < 2)
    matrix2 = matrix1;
end
 
if (nargin < 3)
    if ~exist('n')
        n = mulMatCross(zeros(length(matrix1)),zeros(length(matrix2)));
    end
end

Zmat = ones(2*size(matrix2)-1); % used for a shift right & down - zeroing relavent culmn & row
Zmat(1,:) = NaN;
Zmat(:,1) = NaN;
onesMat = ones(size(matrix2)); % for interior summation
% basically doing this: cov(x,y)/(sigma(x)*sigma(y). see http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
sumXY = circshift(xcorr2(matrix1,matrix2),[1 1]) .* Zmat;
sumX  = circshift(xcorr2(matrix1,onesMat),[1 1]) .* Zmat;
sumY  = circshift(xcorr2(onesMat,matrix2),[1 1]) .* Zmat;
xSq   = circshift(xcorr2(matrix1.^2,onesMat),[1 1]) .* Zmat;
ySq   = circshift(xcorr2(onesMat,matrix2.^2),[1 1]) .* Zmat;
correlationMatrix = (n .* sumXY - sumX .* sumY) ./ ((sqrt(n .* xSq - (sumX).^2) .* sqrt((n+1) .* ySq - (sumY).^2)) + eps);

end

function n_mat = mulMatCross(matrix1,matrix2)
[ma,na] = size(matrix1);
[mb,nb] = size(matrix2);

mc = max([ma+mb-1,ma,mb]);
nc = max([na+nb-1,na,nb]);

m = min(mc,nc);
n_mat = nan(m,m);

i_size = size(matrix2,1);
j_size = size(matrix2,2);

[work_mat,npad_i,npad_j] = pad_edges(matrix1,matrix2,1);

for i = 1:2*size(matrix1,1)-1
    for j = 1:2*size(matrix2,2)-1
        sub_mat = work_mat(npad_i+i-floor(i_size):npad_i+i-1, ...
            npad_j+j-floor(j_size):npad_j+j-1  );
        nan_sub_mat = sub_mat .* matrix2;
        notnan_inds = find(~isnan(nan_sub_mat)); % normalized to the number of nontnan components (average)
        
        n_mat(i,j) = length(notnan_inds);
    end
end

end

function [out_mat,npad_i,npad_j] = pad_edges(mat,h,l)
npad_ij = ceil(size(h)/l);
npad_i = npad_ij(1);
npad_j = npad_ij(2);
in_size = size(mat);
out_size = in_size + [2*npad_i 2*npad_j];
out_mat = nan(out_size);
out_mat(npad_i+1:npad_i+in_size(1),npad_j+1:npad_j+in_size(2)) = mat;

end
