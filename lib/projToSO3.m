function [R] = projToSO3(A)
%PROJTOSO3 Summary of this function goes here
%   Project a matrix A to SO(3) to get a rotation matrix
[U,S,V] = svd(A);
if det(U*V') >= 0
    R = U*V';
else
    R = U*diag([1,1,-1])*V';
end
end

