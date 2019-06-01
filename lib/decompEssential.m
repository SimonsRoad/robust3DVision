function [R,t] = decompEssential(E)
%DECOMPESSENTIAL Summary of this function goes here
%   Decompose Essential matrix to get rotation and translation
[U,S,V]=svd(E);
Rz_p=[cos(pi/2), -sin(pi/2), 0;
    sin(pi/2), cos(pi/2), 0;
    0, 0, 1];
Rz_n=[cos(-pi/2), -sin(-pi/2), 0;
    sin(-pi/2), cos(-pi/2), 0;
    0, 0, 1];
R1=U*Rz_p*V';
R2=U*Rz_n*V';
t1=vmap(U*Rz_p*S*U');
t2=vmap(U*Rz_n*S*U');

R(:,:,1)=R1;
t(:,:,1)=t1;

R(:,:,2)=R1;
t(:,:,2)=t2;

R(:,:,3)=R2;
t(:,:,3)=t1;

R(:,:,4)=R2;
t(:,:,4)=t2;
end

