% % Copyright 2008 Information Systems Lab, ECSE Department, Rensselaer
% % Polytechnic Institute (RPI). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
%     [n,ns,np,ni] = BNET_SIZES(bnet)
%
% function receives a bnet (BNT package) and returns
%   n: number of nodes
%   ns: number of states of each node
%   np: number of distinct parent configurations of each node
%   ni: number of elements in the CPT of each node
function [n,ns,np,ni] = bnet_sizes(bnet)
  n = size(bnet.dag,1);
  ns = bnet.node_sizes(:)';
  for i = 1:n
    % to ensure that node_size of non-discrete nodes is one. This
    % is used in the constrained learning to verify whether a node
    % is discrete or not and to not process non-discrete nodes.
    if ~any(i == bnet.dnodes)
      ns(i)=1;
    end
  end
  % np(i) is the number of parent configurations of node i
  np = ones(1,n);
  for i = 1:n
    for j = 1:n
      if bnet.dag(j,i) == 1
          np(i) = np(i) * ns(j);
      end
    end
  end
  ni = np.*ns;
end
