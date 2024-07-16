function result = meanExclDiag(sqMatrix)
% Get the mean of a square matrix excluding the diagonal, and
% excluding NaN's
dim = size(sqMatrix);
if (length(dim)~=2)||(dim(1)~=dim(2))
    disp(dim);
    error('input is not a square matrix');
end
dim = dim(1);
result = nanmean(sqMatrix(~eye(dim)));
end