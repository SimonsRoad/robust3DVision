% This script solves two view geometry using convex relaxation.
% It first samples groud truth rotation and translations, as well as pairs
% of bearing vectors, then relaxes the problem into an convex optimization
% and solves it.
% Author: Heng Yang           
% Last update: 11/19/2018
%% start clean
clc
clear
close all
cvx_setup
clc

addpath(genpath('./lib'))

NUM_ITERATIONS = 1; % run SDP multiple times
nrPoints = 50;
syntheticOptions.Noise = 0; % in pixel
syntheticOptions.FOV = 100; % in degress
syntheticOptions.MaxParallax = 2; % in meters
syntheticOptions.OutliersRatio=0.0;
syntheticOptions.PlotResiduals = false; % true to turn plot of residuals
syntheticOptions.PlotCameras = false;
view_central_experiment_opt = toVarargin(syntheticOptions);


QCQPOptions.CVXDisp = false;
QCQPOptions.QCQPDisp = false;

finalErrorQCQP = zeros(NUM_ITERATIONS, 3);
finalTightness = zeros(NUM_ITERATIONS, 1);
for itr=1:NUM_ITERATIONS
    %% generates relative rotation, tranlsation and bearing vectors
    [bearing1, bearing2, t_gt, R_gt, outliers_idx] = create2ViewCentralExperiment(nrPoints, view_central_experiment_opt{:});
    %% convex relaxation using QCQP
    [RStarQCQP,tStarQCQP,ZStarQCQP,recoverStatusQCQP,...
        cvxStatusQCQP,pStarRelaxQCQP,tightnessQCQP,rankZ] = ...
        qcqprelaxation(bearing1, bearing2, QCQPOptions);
    
    if recoverStatusQCQP
        rotErrorQCQP = zeros(1,size(RStarQCQP, 2));
        transErrorQCQP = zeros(1,size(RStarQCQP, 2));
        combErrorQCQP = zeros(1,size(RStarQCQP, 2));

        for i=1:size(RStarQCQP, 2)
            rotErrorQCQP(i) = evaluateRotationError(RStarQCQP{i}, R_gt);
            transErrorQCQP(i) = norm(tStarQCQP{i} - t_gt/norm(t_gt));
            combErrorQCQP(i) = sqrt(rotErrorQCQP(i)^2+transErrorQCQP(i)^2);
        end
        [~, bestIdx] = min(combErrorQCQP);
        bestRQCQP = RStarQCQP{bestIdx};
        besttQCQP = tStarQCQP{bestIdx};
        bestErrorQCQP = [rotErrorQCQP(bestIdx), ...
            transErrorQCQP(bestIdx), combErrorQCQP(bestIdx)];
        
        fprintf('QCQP Relaxation tightness: %s, rankZ: %g, rotation error: %g, translation error: %g.\n',...
            string(tightnessQCQP),rankZ,bestErrorQCQP(1),bestErrorQCQP(2));
        finalErrorQCQP(itr,:) = bestErrorQCQP;
        finalTightness(itr,:) = tightnessQCQP;
    end
end





