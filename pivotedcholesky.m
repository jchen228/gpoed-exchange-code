function [F,S,i] = pivotedcholesky(kernel,  pts, k, method)
S = zeros(1,k);
    
    [n,~] = size(pts);
    % diags = kernel(0)*ones(n,1);
    diags = kernel(pts(1,:),pts(1,:))*ones(n,1);
    F = zeros(n, k);
    for i = 1:k
        switch method
            case 'greedy'
                [~,ind] = max(diags);
                si = ind(1);

                % si = i; % no pivot

                S(i) = si;

            case 'rpcholesky'
                % tag = diags./sum(diags);
                % tag(1:10)
                if ~any(diags./sum(diags) > 0)
                    disp(['sensor selection terminated early at k=', num2str(i)])
                    break
                else
                    si = randsample(n, 1, true, diags./sum(diags));
                end
                S(i) = si;
            otherwise
                disp('Method not implemented');
        end
        
        % g = kernel(pdist2(pts, pts(si,:)));
        g = kernel(pts,pts(si,:));
        g = g - F(:,1:(i-1))*F(si,1:(i-1))';
        F(:,i) = g/sqrt(g(si)); % makes col norm of R = 1 when done
        diags = diags - F(:,i).^2;
        % diags(find(diags <0)) = 0;
        diags = diags.*(diags>=0);

        % vecnorm(F')
        % I = eye(n);
        % P = I(:,nonzeros(S));
        % R = F'*P

    end
end