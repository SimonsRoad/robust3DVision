function [R,t,nSteps] = qcqprelaxation_ransac(bearing1,bearing2,...
    nQCQP,R_gt,t_gt,inlierResThresh,maxItr)
%QCQPRELAXATION_RANSAC Use QCQP relaxation inside RANSAC

NUM_POINTS=size(bearing1,2);
if NUM_POINTS < 11
    error('QCQP RANSAC requires at least 11 points')
end

QCQP_QUIET=true;

fullIndex=1:NUM_POINTS;
inlierThresh=0.5;
inlierRatio=0;
itr=0;
while inlierRatio < inlierThresh && itr < maxItr
    randIndex=datasample(fullIndex,nQCQP);
    randBearing1=bearing1(:,randIndex);
    randBearing2=bearing2(:,randIndex);
    % Run QCQP
    [RStarQCQP,tStarQCQP,ZStarQCQP,recoverStatusQCQP,...
        cvxStatusQCQP,pStarRelaxQCQP,tightnessQCQP] = ...
        qcqprelaxation(randBearing1, randBearing2, QCQP_QUIET);
    if isnan(pStarRelaxQCQP)
        disp('WARN: CVX failed, skipping.')
    else
        if ~recoverStatusQCQP
            disp('WARN: QCQP is not tight, skipping')
        else
            rotErrorQCQP = zeros(1,size(RStarQCQP, 2));
            transErrorQCQP = zeros(1,size(RStarQCQP, 2));
            combErrorQCQP = zeros(1,size(RStarQCQP, 2));
            for i=1:size(RStarQCQP, 2)
                rotErrorQCQP(i) = norm(RStarQCQP{i} - R_gt, 'fro');
                transErrorQCQP(i) = norm(tStarQCQP{i} - t_gt/norm(t_gt));
                combErrorQCQP(i) = sqrt(rotErrorQCQP(i)^2+transErrorQCQP(i)^2);
            end
            [~, bestIdx] = min(combErrorQCQP);
            bestRQCQP = RStarQCQP{bestIdx};
            besttQCQP = tStarQCQP{bestIdx};
            % calculate inlier ratio using residual 
            besttQCQP=besttQCQP/norm(besttQCQP);
            bestrQCQP=reshape(bestRQCQP,[9,1]);
            residuals=zeros(1,NUM_POINTS);
            for i=1:NUM_POINTS
                f1=bearing1(:,i);
                f2=bearing2(:,i);
                A=kron(f2',hatmap(f1));
                residuals(i)=(besttQCQP'*A*bestrQCQP)^2;
            end
            inlierRatio=sum(residuals<inlierResThresh)/NUM_POINTS;
            inlierIndex=find(residuals<inlierResThresh);
        end
    end
    itr=itr+1;
end

inlierBearing1=bearing1(:,inlierIndex);
inlierBearing2=bearing2(:,inlierIndex);
[RStarQCQP,tStarQCQP,ZStarQCQP,recoverStatusQCQP,...
        cvxStatusQCQP,pStarRelaxQCQP,tightnessQCQP] = ...
        qcqprelaxation(inlierBearing1, inlierBearing2, QCQP_QUIET);
rotErrorQCQP = zeros(1,size(RStarQCQP, 2));
transErrorQCQP = zeros(1,size(RStarQCQP, 2));
combErrorQCQP = zeros(1,size(RStarQCQP, 2));
for i=1:size(RStarQCQP, 2)
    rotErrorQCQP(i) = norm(RStarQCQP{i} - R_gt, 'fro');
    transErrorQCQP(i) = norm(tStarQCQP{i} - t_gt/norm(t_gt));
    combErrorQCQP(i) = sqrt(rotErrorQCQP(i)^2+transErrorQCQP(i)^2);
end
[~, bestIdx] = min(combErrorQCQP);
bestRQCQP = RStarQCQP{bestIdx};
besttQCQP = tStarQCQP{bestIdx};    

R=bestRQCQP;
t=besttQCQP;
nSteps=itr;
end

