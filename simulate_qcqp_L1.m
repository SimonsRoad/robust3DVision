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

NUM_ITERATIONS = 1; % run SDP multiple times
NOISE_LEVEL = 1; % in pixel
FIELD_OF_VIEW = 100; % in degress
NUM_POINTS = 100;
MAX_PARALLAX = 2; % in meters
OUTLIER_RATIO = 0.0;
DEBUG_SYNTHETIC = false; % true to turn plot on
SDP_QUIET = true; % false to turn cvx output on
QCQP_QUIET = false;

finalErrorQCQP = zeros(NUM_ITERATIONS, 3);
finalTightness = zeros(NUM_ITERATIONS, 1);
for itr=1:NUM_ITERATIONS
    %% generates relative rotation, tranlsation and bearing vectors
    [bearing1, bearing2, t_gt, R_gt] = create2ViewCentralExperiment(...
        NOISE_LEVEL, NUM_POINTS, FIELD_OF_VIEW, MAX_PARALLAX,...
        OUTLIER_RATIO, DEBUG_SYNTHETIC);   
    disp('true R')
    disp(R_gt)
    disp('true t')
    disp(t_gt/norm(t_gt))
    
    
    %% convex relaxation using QCQP
    [RStarQCQP,tStarQCQP,ZStarQCQP,recoverStatusQCQP,...
        cvxStatusQCQP,pStarRelaxQCQP,tightnessQCQP] = ...
        qcqprelaxation_L1(bearing1, bearing2, QCQP_QUIET);
    
    disp(cvxStatusQCQP)
    disp(pStarRelaxQCQP)
    disp('QCQP relax tightness:')
    disp(tightnessQCQP)
    %{
    if recoverStatusQCQP
        rotErrorQCQP = zeros(1,size(RStarQCQP, 2));
        transErrorQCQP = zeros(1,size(RStarQCQP, 2));
        combErrorQCQP = zeros(1,size(RStarQCQP, 2));
        disp('------------recovered R and t----------------')
        for i=1:size(RStarQCQP, 2)
            disp(RStarQCQP{i})
            disp(tStarQCQP{i})
        end
        for i=1:size(RStarQCQP, 2)
            rotErrorQCQP(i) = norm(RStarQCQP{i} - R_gt, 'fro');
            transErrorQCQP(i) = norm(tStarQCQP{i} - t_gt/norm(t_gt));
            combErrorQCQP(i) = sqrt(rotErrorQCQP(i)^2+transErrorQCQP(i)^2);
        end
        [~, bestIdx] = min(combErrorQCQP);
        bestRQCQP = RStarQCQP{bestIdx};
        besttQCQP = tStarQCQP{bestIdx};
        bestErrorQCQP = [rotErrorQCQP(bestIdx), ...
            transErrorQCQP(bestIdx), combErrorQCQP(bestIdx)];
        disp('------------Best R estimation:------------')
        disp(bestRQCQP)
        disp('------------Best t estimation:------------')
        disp(besttQCQP)
        disp('------------Best error:-------------------')
        disp(bestErrorQCQP) 
        disp(itr)
        finalErrorQCQP(itr,:) = bestErrorQCQP;
        finalTightness(itr,:) = tightnessQCQP;
    end 
    %}
    
 
end

% disp('finalErrorQCQP')
% disp(finalErrorQCQP)





