% % Copyright 2008 Information Systems Lab, ECSE Department, Rensselaer
% % Polytechnic Institute (RPI). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% This function has two uses:
%
% FIRST:        p = BNET_PARAM_POS(bn,xi,statexi,varargin)
% xi is the variable number, from 1 to nvar
% statexi is the state of xi
% varargin may have two notations:
% 1) a single number indicating the "binary" representation of the parent
% configuration of xi (use 1 in case of no parents). E.g.
%       bnet_param_pos(bnet,1,1,3)
%    means that x1=1, and parents(x1)=3 (third configuration w.r.t. bnet order)
% 2) a sequence of pair of integer (parent,state) where
% parent must be one of xi parents and state is a valid state for parent.
% The sequence must be preceded by the string 'parents'. E.g.
%       bnet_param_pos(bnet,7,2,'parents',4,1,3,1,6,2)
%    means that x7=2, and parents are x4=1, x3=1 and x6=2.
%    A single vector (or cell vector) of size n (n number of nodes of bnet)
%    is also allowed. In this case, the parent configuration will be
%    extracted from this vector (nodes that are not parents may have any
%    value in the vector, as they will not be used).
%
% % SECOND(reverse mode): (THIS VERSION IS RELEVANT ONLY FOR "INTERNAL" 
% %                        IMPLEMENATION ISSUES)
%       [p,reti,retk] = BNET_PARAM_POS(bn,xi,statexi,varargin)
% % xi is the string 'reverse'
% % if there is only three args, then statexi is the position of a
% % parameter in the vector of parameters for which we want to know
% % the corresponding element in the network. If two extra args are
% % given, then
% % statexi is the number of a variable and the 4th arg is a single
% % number parameter position for which we 
% % want to know the value of parent states of statexi in the network
% % for that given configuration. The return value is p, a cell array
% % with the variable numbers for the parents and their associated
% % states.
% %
% % bn has two different meanings: (i) it is a bnet, (ii) it is a cell array
% % with 5 positions as [bnet, #ofnodes, #statespernode,
% % #parentconfigspernode, sizeCPTpernode]. This strange format is related to
% % performance issue. If we use (i), the function is much slower because it
% % has to evaluate each one of those values.
function [p,reti,retk] = bnet_param_pos(bn,xi,statexi,varargin)
    if iscell(bn)
        bnet = bn{1};
        nn = bn{2};
        ns = bn{3};
        np = bn{4};
        ni = bn{5};
    else
        [nn,ns,np,ni] = bnet_sizes(bn);
        bnet = bn;
    end
    if strcmp(xi,'reverse')
        i = statexi;
        if nargin > 3
            if i < 1 || i > nn
                error('Invalid position to reverse in bnet_param_pos');
            end
            j = varargin{1};
            if j < 1 || j > ni(i)
                error('Invalid position to reverse in bnet_param_pos');
            end
            j = j-1;
        else
            if statexi < 1 || statexi > sum(ni)
                error('Invalid position to reverse in bnet_param_pos');
            end
            sni = 0;
            i = 1;
            while statexi > sni+ni(i)
                sni=sni+ni(i);
                i = i + 1;
            end
            j = statexi-sni-1;
            reti = mod(j,np(i))+1;
            retk = div(j,np(i))+1;
        end
        j = mod(j,np(i));
        reti = i;
%        p = cell(1,max(find(bnet.dag(:,reti))));
        p = cell(1,nn);
        for i = 1:nn
            if bnet.dag(i,reti)
                p{i} = mod(j,ns(i))+1;
                j = div(j,ns(i));
            end
        end
    else
        reti=0;
        retk=0;
        if xi < 1 || xi > nn || statexi < 1 || statexi > ns(xi)
            fprintf('x%d state %d\n', xi, statexi);
            error('Invalid (2nd or 3rd) parameter to bnet_param_pos');
        end
        parentconfig = 1;
        if nargin > 3
            if strcmp(varargin{1},'parents')
                if numel(varargin{2})>1 || iscell(varargin{2})
                    parents = varargin{2};
                    if iscell(parents)
                        temp = parents;
                        parents = zeros(1,nn);
                        for i = 1:numel(temp)
                            if temp{i}>0
                                parents(i)=temp{i};
                            end
                        end
                    end
                else
                    n = 2;
                    parents = zeros(nn,1);
                    while n < nargin-3
                        pxi = varargin{n};
                        statepxi = varargin{n+1};
                        if pxi < 1 || pxi > nn || ~bnet.dag(pxi,xi)
                            fprintf('x%d is not parent of x%d. Ignoring...\n', pxi, xi);
                            %                    error('Invalid parameter to bnet_param_pos');
                        end
                        if statepxi < 1 || statepxi > ns(pxi)
                            fprintf('Invalid state %d for x%d\n', statepxi, pxi);
                            error('Invalid parameter to bnet_param_pos');
                        end
                        parents(pxi)=statepxi;
                        n = n + 2;
                    end
                end
                ap = 1;
                for j = find(bnet.dag(:,xi))'
                    %            for j = 1:nn
                    %                if bnet.dag(j,xi)
                    if parents(j)<1
                        fprintf('State of parent x%d of x%d was not defined\n', j, xi);
                        error('Invalid parameter to bnet_param_pos');
                    end
                    parentconfig = parentconfig + ap*(parents(j)-1);
                    ap = ap * ns(j);
                    %                end
                end
            else
                if varargin{1} < 1 || varargin{1} > np(xi)
                    error('Invalid 4th parameter to bnet_param_pos');
                end
                parentconfig = varargin{1};
            end
        else
            if np(xi)>1
                error('Parents were not defined in bnet_param_pos');
            end
        end
        p = sum(ni(1:xi-1)) + parentconfig + (statexi-1)*np(xi);
    end
end
