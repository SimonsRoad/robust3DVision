function [ stable_rank ] = stableRank( A )
%STABLERANK Calculates numerical/stable rank of matrix
stable_rank = norm(A,'fro')^2 / norm(A,2)^2;
end

