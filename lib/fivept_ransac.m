function [R_est,t_est,nSteps] = fivept_ransac(bearing1,bearing2,...
    R_gt,t_gt,inlierResThresh,maxItr)
%FIVEPT_RANSAC Call opengv five point inside RANSAC

NUM_POINTS=size(bearing1,2);
if NUM_POINTS < 10
    error('5-POINT RANSAC requires at least 10 points')
end

fullIndex=1:NUM_POINTS;
inlierThresh=0.5;
inlierRatio=0;
itr=0;
nPoints=5;

t_gt=t_gt/norm(t_gt);
while inlierRatio < inlierThresh && itr < maxItr
    randIndex=datasample(fullIndex,nPoints);
    % call opengv five point method
    X_five_pt=opengv('fivept_nister',randIndex,bearing1,bearing2);
    if size(X_five_pt,3)==0
        disp('Warning: opengv 5-point returns empty solution, skipping.')
    else
        % transform essential to rotation matrices
        R_five_pt=transformEssentials(X_five_pt);
        % take the best R compared to gt
        temp=zeros(1,size(R_five_pt,3));
        for k=1:size(R_five_pt,3)
            temp(k)=evaluateRotationError(R_gt,R_five_pt(:,:,k));
        end
        [~,minIdx]=min(temp);
        R_five_pt_best=R_five_pt(:,:,minIdx);
        E_five_pt_best=X_five_pt(:,:,ceil(minIdx/2));
        % use best R to get t
        t_five_pt_best=vmap(E_five_pt_best*R_five_pt_best');
        t_five_pt_best=t_five_pt_best/norm(t_five_pt_best);
        if norm(-t_five_pt_best-t_gt) < norm(t_five_pt_best-t_gt)
            t_five_pt_best=-t_five_pt_best;
        end
        % use best R and t to calculate residuals and # of inliers
        r_five_pt_best=reshape(R_five_pt_best,[9,1]);
        residuals=zeros(1,NUM_POINTS);
        for pt=1:NUM_POINTS
            f1=bearing1(:,pt);
            f2=bearing2(:,pt);
            A=kron(f2',hatmap(f1));
            residuals(pt)=(t_five_pt_best'*A*r_five_pt_best)^2;
        end
        inlierRatio=sum(residuals<inlierResThresh)/NUM_POINTS;
        inlierIndex=find(residuals<inlierResThresh);
    end
    itr=itr+1; 
end
nSteps=itr;
if nSteps > (maxItr-5)
    disp('WARN: FIVE POINT RANSAC did not converge within maximum iterations')
end
% use nonlinear optimization to refine result
X_5pt_nonlin_opt=opengv('rel_nonlin_central',inlierIndex,...
            bearing1,bearing2,[R_five_pt_best,t_five_pt_best]);
R_est=X_5pt_nonlin_opt(:,1:3);
t_est=X_5pt_nonlin_opt(:,4)/norm(X_5pt_nonlin_opt(:,4));

end

