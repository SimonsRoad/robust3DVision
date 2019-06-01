% This script solves two view geometry using convex relaxation.
% It first samples groud truth rotation and translations, as well as pairs
% of bearing vectors, then relaxes the problem into an convex optimization
% and solves it.
% Author: Heng Yang           
% Last update: 11/26/2018
% simulate QCQP with different outlier ratios, 100 observations, 30
% iterations
%% start clean
clc
clear
close all

NUM_ITERATIONS = 30; % run SDP multiple times
NOISE_LEVEL = 0:2:20; % in pixel
FIELD_OF_VIEW = 100; % in degress
NUM_POINTS = 50;
MAX_PARALLAX = 2; % in meters
OUTLIER_RATIO = 0.0;
DEBUG_SYNTHETIC = false; % true to turn plot on
SDP_QUIET = true; % false to turn cvx output on
QCQP_QUIET = true;

performance = {};
finalErrorQCQP = zeros(NUM_ITERATIONS, 3);
finalTightness = zeros(NUM_ITERATIONS, 1);
for noise_itr=1:length(NOISE_LEVEL)
    noise_level = NOISE_LEVEL(noise_itr);
    for itr=1:NUM_ITERATIONS
        %% generates relative rotation, tranlsation and bearing vectors
        [bearing1, bearing2, t_gt, R_gt] = create2ViewCentralExperiment(...
            noise_level, NUM_POINTS, FIELD_OF_VIEW, MAX_PARALLAX,...
            OUTLIER_RATIO, DEBUG_SYNTHETIC);   
        disp('true R')
        disp(R_gt)
        disp('true t')
        disp(t_gt/norm(t_gt))


        %% convex relaxation using QCQP
        [RStarQCQP,tStarQCQP,ZStarQCQP,recoverStatusQCQP,...
            cvxStatusQCQP,pStarRelaxQCQP,tightnessQCQP] = ...
            qcqprelaxation(bearing1, bearing2, QCQP_QUIET);

        disp(cvxStatusQCQP)
        disp(pStarRelaxQCQP)
        disp('QCQP relax tightness:')
        disp(tightnessQCQP)
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
    end
    performance{noise_itr,1} = finalErrorQCQP;
    performance{noise_itr,2} = finalTightness;
end

for i=1:length(NOISE_LEVEL)
pose_error(:,i) = performance{i,1}(:,1);
end
boxplot(pose_error)

