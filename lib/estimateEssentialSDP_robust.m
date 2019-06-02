function [R_est,t_est,E_est,f_sdp,f_round,Z_sdp,thetas_sdp,solverOutput] = estimateEssentialSDP_robust(bearing1,bearing2,varargin)
% Robust version of estimateEssentialSDP, using Truncated Least Squares
% cost function, with binary cloning and SDP relaxation
% Author: Heng Yang
% Last update: 06/02/2019
% Massachusetts Institute of Technology

