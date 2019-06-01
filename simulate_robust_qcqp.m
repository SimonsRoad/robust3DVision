% simulate robust QCQP using Geman-McClure function
% solve using graduated non-convexity and QCQP as base solver
%
% Author: Heng Yang
% Massachusetts Institute of Technology
% Feb. 02, 2019

clc
clear
close all
cvx_setup
clc

addpath(genpath('./lib'))

% options for generating synthetic data
nrPoints = 50;
syntheticOptions.Noise = 0; % in pixel
syntheticOptions.FOV = 100; % in degress
syntheticOptions.MaxParallax = 2; % in meters
% syntheticOptions.outlierRatio = 1/syntheticOptions.nrPoints;
syntheticOptions.OutliersRatio=0.2;
syntheticOptions.PlotResiduals = false; % true to turn plot of residuals
syntheticOptions.PlotCameras = false;
view_central_experiment_opt = toVarargin(syntheticOptions);

% solve options for least squares qcqp
QCQPOptions.CVXDisp=false;
QCQPOptions.QCQPDisp=false;

% solve options for robust approach
robustQCQPOptions.CVXDisp = false;
robustQCQPOptions.QCQPDisp = false;
robustQCQPOptions.robustDisp = true;
robustQCQPOptions.maxSteps=1e4;
robustQCQPOptions.initialMu=10;
robustQCQPOptions.divFactor=0.9;
robustQCQPOptions.stopThresh=1e-4;

%% create synthetic data
[bearing1, bearing2, t_gt, R_gt, outlier_idx] = create2ViewCentralExperiment(nrPoints, view_central_experiment_opt{:});

t_gt_norm=norm(t_gt);
t_gt=t_gt/t_gt_norm;

%% solve using least squares QCQP
fprintf('===================== Least Squares QCQP ==================\n')
[RStarQCQP,tStarQCQP,ZStarQCQP,recoverStatusQCQP,...
        cvxStatusQCQP,pStarRelaxQCQP,tightnessQCQP,rankZ] = ...
        qcqprelaxation(bearing1, bearing2, QCQPOptions);
for i=1:size(RStarQCQP, 2)
    rotErrorQCQP(i) = evaluateRotationError(RStarQCQP{i}, R_gt);
    transErrorQCQP(i) = norm(tStarQCQP{i} - t_gt/norm(t_gt));
end 
fprintf('QCQP rotation error: '); fprintf('%g ',rotErrorQCQP); fprintf('\n');
fprintf('QCQP translation error: '); fprintf('%g ',transErrorQCQP); fprintf('\n'); 

%% solve using least squares QCQP with ground truth outlier idices
fprintf('===================== Least Squares QCQP with outliers removed ==================\n')
weights_gt=ones(1,nrPoints);
weights_gt(outlier_idx)=0;
bearing1_gt=weights_gt.*bearing1;
bearing2_gt=weights_gt.*bearing2;
[RStarInlier,tStarInlier,ZStarInlier,recoverStatusInlier,...
        cvxStatusInlier,pStarRelaxInlier,tightnessInlier,rankZInlier] = ...
        qcqprelaxation(bearing1_gt, bearing2_gt, QCQPOptions);
for i=1:size(RStarInlier, 2)
    rotErrorInlier(i) = evaluateRotationError(RStarInlier{i}, R_gt);
    transErrorInlier(i) = norm(tStarInlier{i} - t_gt/norm(t_gt));
end
essentialInlier = hatmap(tStarInlier{1}) * RStarInlier{1};
costInlier=0;
for i=1:nrPoints
	f1=bearing1(:,i);
	f2=bearing2(:,i);
	costInlier=costInlier+weights_gt(i)*(f1'*essentialInlier*f2)^2;
end
fprintf('Inlier QCQP rotation error: '); fprintf('%g ',rotErrorInlier); fprintf('\n');
fprintf('Inlier QCQP translation error: '); fprintf('%g ',transErrorInlier); fprintf('\n'); 
fprintf('Inlier QCQP cost=%g.\n',costInlier);

%% solve using robust approach with Geman-McClure cost function
fprintf('===================== Robust QCQP ==================\n')
[RStarRobust,tStarRobust,EStarRobust,nrItrsRobust,weightsRobust,muRobust,robustQCQPOptions] = ...
    robustQCQP(bearing1,bearing2,robustQCQPOptions);
for i=1:size(RStarRobust, 2)
    rotErrorRobust(i) = evaluateRotationError(RStarRobust{i}, R_gt);
    transErrorRobust(i) = norm(tStarRobust{i} - t_gt/norm(t_gt));
end    

fprintf('Robust approach rotation error: '); fprintf('%g ',rotErrorRobust); fprintf('\n');
fprintf('Robust approach translation error: '); fprintf('%g ',transErrorRobust); fprintf('\n');
figure;
bar(weightsRobust);




