function [result,cvxStatus,pStarRelax] = sdprelaxation(bearing1,bearing2,cvxQuiet)
%SDPRELAXATION Summary of this function goes here
%   Relax two view geometry into SDP
%   Need redundant constraints, otherwise stuck in trivial solution
nSamples = size(bearing1, 2);
if cvxQuiet
    disp('set quiet')
    cvx_begin quiet
end
cvx_begin sdp
cvx_precision best
variable X(12, 12) symmetric
expression temp(3,3)
expression y(nSamples)
for i=1:nSamples
    f2 = bearing2(:, i);
    f1 = bearing1(:, i);
    F1 = hatmap(f1);
    for j = 1:3
        temp(j,:) = F1(j,:) * X((1+3j):(3+3j), 1:3);
    end
    y(i) = abs(trace(f2 * sum(temp, 1)));
end
%     minimize(sdpcost(bearing1, bearing2, X))
minimize(sum(y))
subject to
X >= 0
X(1:3, 1:3) == eye(3)
trace(X(4:end, 4:end)) == 3
cvx_end

result = X;
cvxStatus = cvx_status;
pStarRelax = cvx_optval;

end

