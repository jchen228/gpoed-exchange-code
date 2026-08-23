function [p, ld_exchange] = dist_sel5(x, k, bucket_schedule, sig_n, sig_f, ls, f, kernel_fn)
% DIST_SEL5 Multi-round hierarchical sensor selection with lineage-preserving bucketing
%
% 
%   - Round 1 assigns sensors to buckets via i.i.d. uniform draw.
%   - Rounds 2. never re-cluster. They merge GROUPS of the previous
%     round's buckets (bucket_schedule(r-1)/bucket_schedule(r) consecutive
%     buckets per new bucket) and re-optimize each merged group locally,
%     preserving bucket lineage instead of flattening into one pool.
%   - Every round -- intermediate buckets and the final global exchange --
%     uses EXCHANGE_SENSORS_SEQUENTIAL_INCREMENTAL, which maintains the
%     candidate solves incrementally instead of a fresh O(c*k^2) solve
%     every iteration.
%
% bucket_schedule: vector of bucket counts per round, e.g. [64,32,8,4,2,1]
%   Consecutive entries must divide evenly (bucket_schedule(r) must be an
%   integer multiple of bucket_schedule(r+1)) except when transitioning
%   into a terminal round (bucket_schedule(r+1) == 1), which is always
%   allowed regardless of group size.
% kernel_fn: optional, defaults to @gaussKern (isotropic squared
% exponential). Pass @maternKern32/@maternKern52 for ARD Matern -- in
% that case ls should be a 1xD vector rather than a scalar.
if nargin < 8 || isempty(kernel_fn)
    kernel_fn = @gaussKern;
end

K_fun = @(x) kernel_fn(x,sig_f,ls);

n = size(x, 1);
num_rounds = length(bucket_schedule);

for r = 1:num_rounds-1
    if bucket_schedule(r+1) ~= 1 && mod(bucket_schedule(r), bucket_schedule(r+1)) ~= 0
        error('dist_sel5:invalidSchedule', ...
            'bucket_schedule(%d)=%d must be an integer multiple of bucket_schedule(%d)=%d.', ...
            r, bucket_schedule(r), r+1, bucket_schedule(r+1));
    end
end

prev_bucket_results = {};

for r = 1:num_rounds
    num_buckets = bucket_schedule(r);

    if r == 1
        pool_size = n;
    else
        pool_size = sum(cellfun(@numel, prev_bucket_results));
    end

    fprintf('Round %d/%d: %d sensors → %d buckets\n', ...
        r, num_rounds, pool_size, num_buckets);

    %% --- Final round: global exchange ---
    if num_buckets == 1
        if r == 1
            pool_idx = (1:n)';
        else
            pool_idx = cat(1, prev_bucket_results{:});
        end
        pool_size = length(pool_idx);
        x_pool = x(pool_idx, :);

        if pool_size <= k
            p_init = 1:pool_size;
        else
            p_init = round(linspace(1, pool_size, k));
        end

        p_local = exchange_sensors_sequential_incremental(...
            p_init, x_pool, sig_n, sig_f, ls, f, kernel_fn);
        ld_exchange = slogdet(K_fun(x_pool(p_local,:)), sig_n)
        p = pool_idx(p_local);
        return;
    end

    bucket_results = cell(num_buckets, 1);

    if r == 1
        %% --- Round 1: fresh i.i.d. random bucketing ---
        bucket_labels = randi(num_buckets, n, 1);

        parfor b = 1:num_buckets
            idx_bucket = find(bucket_labels == b);
            x_bucket = x(idx_bucket, :);
            bucket_size = length(idx_bucket);

            if bucket_size == 0
                continue;
            end

            k_local = min(k, bucket_size);
            p_init_local = round(linspace(1, bucket_size, k_local));

            p_opt_local = exchange_sensors_sequential_incremental(...
                p_init_local, x_bucket, sig_n, sig_f, ls, f, kernel_fn);

            bucket_results{b} = idx_bucket(p_opt_local(:));
        end
    else
        %% --- Rounds 2+: merge groups of the previous round's buckets ---
        group_size = bucket_schedule(r-1) / num_buckets;

        parfor b = 1:num_buckets
            group_idx = (b-1)*group_size+1 : b*group_size;
            idx_bucket = cat(1, prev_bucket_results{group_idx});
            x_bucket = x(idx_bucket, :);
            bucket_size = length(idx_bucket);

            if bucket_size == 0
                continue;
            end

            k_local = min(k, bucket_size);
            p_init_local = round(linspace(1, bucket_size, k_local));

            p_opt_local = exchange_sensors_sequential_incremental(...
                p_init_local, x_bucket, sig_n, sig_f, ls, f, kernel_fn);

            bucket_results{b} = idx_bucket(p_opt_local(:));
        end
    end

    prev_bucket_results = bucket_results;

    fprintf('Round %d complete: pool reduced to %d sensors\n', ...
        r, sum(cellfun(@numel, bucket_results)));
end

%% --- Safety Fallback Global Exchange ---
% Runs automatically if the bucket_schedule array did not explicitly terminate in 1
pool_idx = cat(1, prev_bucket_results{:});
pool_size = length(pool_idx);
x_pool = x(pool_idx, :);

if pool_size <= k
    p_init = 1:pool_size;
else
    p_init = round(linspace(1, pool_size, k));
end

p_local = exchange_sensors_sequential_incremental(p_init, x_pool, sig_n, sig_f, ls, f, kernel_fn);
ld_exchange = slogdet(K_fun(x_pool(p_local,:)), sig_n)
p = pool_idx(p_local);

end
