% % RP Cholesky function - with option for using tolerance or specified K
% % INPUT:
% % % number of "sensors" N
% % % grid information x: N x D
% % % sigma variance
% % % sig_n noise
% % % varargin: sets whether we use a tolerance based on trace norm or on a
% % fixed number of sensors
% 
% % OUTPUT:
% % % selection operator S (pivot set)
% % % matrix F: n x k defining Nystrom spproximation A_hat = FF'


function [S, F] = rpchol_tol(Kf,varargin)
% ls = 1; sig_n = 0; sig_f = 1;

    if isa(Kf,'funMat')
        N = max(Kf.Asize);
        
        diag_KXX_hat = zeros(N,1);
        F = []; d = diagkern(Kf)'; % diagonal elements of KXX
        d1 = d; % save diagonal elements of kernel matrix 
        prob = d/sum(d);
        S = []; % initialize
        
        % i = 1 step
        i = 1;
        S(i) = datasample(1:N,1,'Replace',false,'Weights',prob);
        % g = kxx_col(S(1), x, ls, sig_n); % no removal step
        % g = Kf*(double(1:N == S(i))');
        g = rowcolkern(Kf,1:N,S(i));
        F(:,1) = g/sqrt(g(S(1)));
        d = d - abs(F(:,1)).^2;
        d = d.*(d>=0); % ensure d is non-negative
        
        diag_KXX_hat = sum(F.^2,2); 
        
        if varargin{1} < 1
            disp('epsilon tolerance')
            epsilon = varargin{1};
        
            while  sum(d1 - diag_KXX_hat) > epsilon*sum(d1) 
                
                i = i+1
                prob = d/sum(d); % rewrite probability
                S(i) = datasample(1:N,1,'Replace',false,'Weights',prob);
                % g = kxx_col(S(i), x, ls, sig_n,sig_f); % find s-th column of A
                % g = Kf*(double(1:N == S(i))');
                g = rowcolkern(Kf,1:N,S(i));
                g = g - F(:,1:(i-1))*F(S(i),1:(i-1))'; % remove overlap with previously chosen columns
                F(:,i) = g/sqrt(g(S(i))); % update approximation
                d = d - abs(F(:,i)).^2;
                d = d.*(d>=0);
            
                diag_KXX_hat = sum(F.^2,2); 
            end
            sum(d1 - diag_KXX_hat)/sum(d1) 
        
        else
            m = varargin{1};
            disp('set number of sensors')
        
            for i = 2:m
                prob = d/sum(d); % rewrite probability
                S(i) = datasample(1:N,1,'Replace',false,'Weights',prob);
                % g = kxx_col(S(i), x, ls, sig_n,sig_f); % find s-th column of A
                % g = Kf*(double(1:N == S(i))');
                g = rowcolkern(Kf,1:N,S(i));
                g = g - F(:,1:(i-1))*F(S(i),1:(i-1))'; % remove overlap with previously chosen columns
                F(:,i) = g/sqrt(g(S(i))); % update approximation
                d = d - abs(F(:,i)).^2;
                d = d.*(d>=0);
            end
        
        end

    elseif isa(Kf, "function_handle")
        disp('function handle')
        x = varargin{2};
        N = length(x);
                
        diag_KXX_hat = zeros(N,1);
        F = []; d = Kf(x,'diag'); % diagonal elements of KXX
        d1 = d; % save diagonal elements of kernel matrix 
        prob = d/sum(d);
        S = []; % initialize
        
        % i = 1 step
        i = 1;
        S(i) = datasample(1:N,1,'Replace',false,'Weights',prob);
        % g = kxx_col(S(1), x, ls, sig_n); % no removal step
        % g = Kf(1:N,S(i));
        g = Kf(x,x(S(i),:));
        F(:,1) = g/sqrt(g(S(1)));
        d = d - abs(F(:,1)).^2;
        d = d.*(d>=0); % ensure d is non-negative
        
        diag_KXX_hat = sum(F.^2,2); 
        
        if varargin{1} < 1
            disp('epsilon tolerance')
            epsilon = varargin{1};
        
            while  sum(d1 - diag_KXX_hat) > epsilon*sum(d1) 
                
                i = i+1
                prob = d/sum(d); % rewrite probability
                S(i) = datasample(1:N,1,'Replace',false,'Weights',prob);
                % g = kxx_col(S(i), x, ls, sig_n,sig_f); % find s-th column of A
                % g = Kf(1:N,S(i));
                g = Kf(x,x(S(i),:));
                g = g - F(:,1:(i-1))*F(S(i),1:(i-1))'; % remove overlap with previously chosen columns
                F(:,i) = g/sqrt(g(S(i))); % update approximation
                d = d - abs(F(:,i)).^2;
                d = d.*(d>=0);
            
                diag_KXX_hat = sum(F.^2,2); 
            end
            sum(d1 - diag_KXX_hat)/sum(d1) 
        
        else
            m = varargin{1};
            disp('set number of sensors')
        
            for i = 2:m
                prob = d/sum(d); % rewrite probability
                S(i) = datasample(1:N,1,'Replace',false,'Weights',prob);
                % g = kxx_col(S(i), x, ls, sig_n,sig_f); % find s-th column of A
                % g = Kf(1:N,S(i));
                g = Kf(x,x(S(i),:));
                g = g - F(:,1:(i-1))*F(S(i),1:(i-1))'; % remove overlap with previously chosen columns
                F(:,i) = g/sqrt(g(S(i))); % update approximation
                d = d - abs(F(:,i)).^2;
                d = d.*(d>=0);
            end
        
        end
    else
        N = length(Kf);
                
        diag_KXX_hat = zeros(N,1);
        F = []; d = diag(Kf); % diagonal elements of KXX
        d1 = d; % save diagonal elements of kernel matrix 
        prob = d/sum(d);
        S = []; % initialize
        
        % i = 1 step
        i = 1;
        S(i) = datasample(1:N,1,'Replace',false,'Weights',prob);
        % g = kxx_col(S(1), x, ls, sig_n); % no removal step
        g = Kf(1:N,S(i));
        F(:,1) = g/sqrt(g(S(1)));
        d = d - abs(F(:,1)).^2;
        d = d.*(d>=0); % ensure d is non-negative
        
        diag_KXX_hat = sum(F.^2,2); 
        
        if varargin{1} < 1
            disp('epsilon tolerance')
            epsilon = varargin{1};
        
            while  sum(d1 - diag_KXX_hat) > epsilon*sum(d1) 
                
                i = i+1
                prob = d/sum(d); % rewrite probability
                S(i) = datasample(1:N,1,'Replace',false,'Weights',prob);
                % g = kxx_col(S(i), x, ls, sig_n,sig_f); % find s-th column of A
                g = Kf(1:N,S(i));
                g = g - F(:,1:(i-1))*F(S(i),1:(i-1))'; % remove overlap with previously chosen columns
                F(:,i) = g/sqrt(g(S(i))); % update approximation
                d = d - abs(F(:,i)).^2;
                d = d.*(d>=0);
            
                diag_KXX_hat = sum(F.^2,2); 
            end
            sum(d1 - diag_KXX_hat)/sum(d1) 
        
        else
            m = varargin{1};
            disp('set number of sensors')
        
            for i = 2:m
                prob = d/sum(d); % rewrite probability
                S(i) = datasample(1:N,1,'Replace',false,'Weights',prob);
                % g = kxx_col(S(i), x, ls, sig_n,sig_f); % find s-th column of A
                g = Kf(1:N,S(i));
                g = g - F(:,1:(i-1))*F(S(i),1:(i-1))'; % remove overlap with previously chosen columns
                F(:,i) = g/sqrt(g(S(i))); % update approximation
                d = d - abs(F(:,i)).^2;
                d = d.*(d>=0);
            end
        
        end
    end