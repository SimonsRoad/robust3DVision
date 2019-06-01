function [gam,vars,opt,prog] = sosrelaxation(bearing1,bearing2)
%UNTITLED Summary of this function goes here
%   solve two-view geometry using SOS
nSamples = size(bearing1, 2);
r = sym('r', [9,1], 'real');
R = reshape(r, [3,3]);
t = sym('t', [3,1], 'real');
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


%{
temp = sym('temp', [3,3], 'real');
y = sym('y', [1, nSamples], 'real');

for i=1:nSamples
    f2 = bearing2(:, i);
    f1 = bearing1(:, i);
    F1 = hatmap(f1);
    for j = 1:3
        disp(F1(j,:))
        temp(j,:) = t(j) * F1(j,:)  * R;
%         disp('------temp------')
%         disp(temp(j,:))
    end
    y(i) = trace(f2 * sum(temp, 1))^2; % quadraric version
%     y(i) = trace(f2 * sum(temp, 1)); % linear version
end
sym cost;
cost = sum(y);
disp(symvar(cost))
vars = sym(symvar(cost));
disp(vars)
disp('----------y------------')
disp(y)
disp('---------cost----------')
disp(cost)
%}


% define SOSP
equalities = sym('equalities', [1, 1]);
count = 1;
for i=1:3
    for j=i:3
        if j==i  % Unit norm
            equalities(count) = R(:,i)'*R(:,i) - 1;
            count = count + 1;
        else % othogonal
            equalities(count) = R(:,i)'*R(:,j);
            count = count + 1;
        end
    end
end
equalities(count) = t'*t - 1;

degree = 4;
options = struct('solver', 'sdpt3');

[gam,vars,opt,prog] = findboundSOS(cost,[],...
    equalities,degree,options);
end

