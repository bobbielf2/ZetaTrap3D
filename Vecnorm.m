function vn = vecnorm(A)
% Alternative vecnorm function (in case MATLAB version earlier than R2017b)
% return column-wise 2-norm for the matrix A.
vn = sqrt(sum(A.^2));
