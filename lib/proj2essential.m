function [E] = proj2essential(M)
% project a matrix M to essential manifold

[U,S,V] = svd(M);
sigmas = diag(S);
temp = (sigmas(1) + sigmas(2))/2;
Sproj = diag([temp,temp,0]);

E = U*Sproj*V';
