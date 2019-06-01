function [sol,u,Q] = sosrelaxationYALMIP(bearing1,bearing2,R_gt,t_gt)
%UNTITLED Summary of this function goes here
%   solve two-view geometry using SOS
nSamples = size(bearing1, 2);
r = sdpvar(9,1);
t = sdpvar(3,1);
varVec = [r;t];
R = reshape(r, [3,3]);
X = r * t';
x = reshape(X, [numel(X), 1]);
C = {3,3};
for i=1:3
    for j=1:3
        e_i = zeros(3,1); e_i(i) = 1;
        e_j = zeros(3,1); e_j(j) = 1;
        C_ij = zeros(9,9);
        for k=1:nSamples
            f1 = bearing1(:, k);
            f2 = bearing2(:, k);
            C_ij = C_ij + ...
                kron(f2, cross(f1, e_j))*(kron(f2, cross(f1, e_i)))';
        end
        C{i,j} = C_ij;
    end
end
C_extend = zeros(numel(X), numel(X));
for i=1:3
    for j=1:3
        e_ij = zeros(3,3); e_ij(i,j) = 1;
        C_extend = C_extend + kron(e_ij, C{i,j});
    end
end

cost = x' * C_extend * x;
disp('---------cost----------')
disp(cost)

sdpvar lowerbound

% define SOSP
count = 1;
for i=1:3
    for j=i:3
        if j==i  % Unit norm
            equalities(count,1) = R(:,i)'*R(:,i) - 1;
            count = count + 1;
        else % othogonal
            equalities(count,1) = R(:,i)'*R(:,j);
            count = count + 1;
        end
    end
end
equalities(count,1) = t'*t - 1;

for i=1:length(equalities)
    equalities(end+1,1) = -1*equalities(i,1);
end
degree = 2;
for i=1:length(equalities)
    [multi_s(i,1), multi_c{i}] = polynomial(varVec',degree);
end

F = [ sos(cost-lowerbound-multi_s'*equalities) ];
for i=1:length(equalities)
    F = [F, sos(multi_s(i,1))];
end
disp(F)
decvar = multi_c{1};
for i=2:length(multi_c)
    decvar = vertcat(decvar, multi_c{i});
end
[sol,u,Q] = solvesos(F, -lowerbound, [], [decvar; lowerbound]);
disp('---------Lower bound------------')
value(lowerbound)

r_gt = reshape(R_gt, [9,1]);
X_gt = r_gt * t_gt';
x_gt = reshape(X_gt, [numel(X_gt), 1]);
cost_gt = x_gt' * C_extend * x_gt;
disp('---------ground truth cost----------')
disp(cost_gt)
end

