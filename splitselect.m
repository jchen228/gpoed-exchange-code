function [pk_p] = splitselect(x,p,k,K)
% splitselection function allows you to split the kernel into submatrices
% and then select (for now) through cssp algorithm
% (not made for parallel processing yet)

% NOTE: need to improve partitioining algorithm. 

% how to use 
% % pk_test = splitselect(x,p,k,@(x,x2) K_fun_offdiag(x,x2));

arguments (Input)
    x % spatial data
    p % number of processors
    k % number of available sensors
    K % kernel function
end

arguments (Output)
    pk_p
end

pk_p = [];

% Calculate the size of each submatrix
n = length(x);
group_size = ceil(n/p);
group_lengths = ones(1, p) * group_size;
k_p = floor(k/p);
k_group_lengths = ones(1, p) * k_p;
k_group_lengths(1:mod(k,p)) = k_group_lengths(end) + 1;
% Adjust for remainders
group_lengths(end) = group_lengths(end) + n - sum(group_lengths);
% Partition spatial domain and k available sensors
parts = mat2cell([1:n]', group_lengths);


for i = 1:p
    K_i = K(x(parts{i},:),x(parts{i},:));
    partselection = gks(K_i,k_group_lengths(i)); % partitioned selection 
    pk_p = [pk_p parts{i}(partselection)'];
end

end