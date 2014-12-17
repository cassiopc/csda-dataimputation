% % Copyright 2008 Information Systems Lab, ECSE Department, Rensselaer
% % Polytechnic Institute (RPI), Troy, NY. All rights reserved.
% % http://www.ecse.rpi.edu/homepages/cvrl/
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
%        ok = BNET_CHECK(bnet)
%
% test if probability values in the bnet are correct. They must be between
% zero and one and must sum one. Return 1 if ok and zero otherwise. In that
% case, it writes some messages to explain why it is not ok.
function ok = bnet_check(bnet,eps,verb)

reqfuncs = { 'bnet_sizes', 'bnet_get_params_as_vector', 'bnet_param_pos' };
for i=1:length(reqfuncs)
   if exist(reqfuncs{i})==0
      fprintf('Required function %s not found\n',reqfuncs{i});
      return
   end
end
    [nn,ns,np,ni] = bnet_sizes(bnet);

    for i = 1:nn
      for j = 1:i
        if bnet.dag(i,j) > 0
            fprintf('WARNING: node order does not respect topological order\n');
        end
      end
    end
    if nargin < 2
    % epsilon is good for float point comparisons
        eps = 0.00001;
    end
    if nargin < 3
      verb = 1;
    end
    x = bnet_get_params_as_vector(bnet);
    %disp(['x=' num2str(x)]);
    ok = 1;
    % for each node
    for i = 1:nn
        % for each parent configuration
        for j = 1:np(i)
            s = 0;
            % for each state
            for k = 1:ns(i)
                v = x(bnet_param_pos(bnet,i,k,j));
                % probability value must lie in [0,1]
                if v < -eps || v > 1+eps
                  if verb > 0
                    fprintf(['Check failed: node %d, parent config %d, state %d, value = %.8f is not in [0,1]\n'],i,j,k,v);
                  end
                  ok = 0;
                end
                s = s + v;
            end
            % sum for all states must be one
            if abs(s-1)>eps
              if verb > 0
                fprintf('Check failed: node %d, parent config %d, sum = %.8f\n',i,j,s);
              end
              ok = 0;
            end
        end
    end
end
