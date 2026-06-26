function [max_beta, best_local_idx] = evaluate_block(L_fast, rhs_chunk, alpha_chunk, beta_orig, f)
% EVALUATE_BLOCK Worker routine executing Level 3 BLAS math on a candidate chunk

    % 1. Level 3 BLAS Triangular Solve (Screamingly fast on short & fat matrices)
    b_chunk = L_fast \ rhs_chunk;
    
    % 2. Vectorized element-wise calculations 
    beta_chunk = sqrt(alpha_chunk' - sum(b_chunk.^2, 1));
    
    % 3. Extract the local maximum
    [max_beta, best_local_idx] = max(beta_chunk);
    
    % Note: We avoid returning the full b_chunk matrix to minimize 
    % communication overhead between worker and master.
end