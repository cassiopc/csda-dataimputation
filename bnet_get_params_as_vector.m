% % Copyright 2008 Information Systems Lab, ECSE Department, Rensselaer
% % Polytechnic Institute (RPI). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
%     x = BNET_GET_PARAMS_AS_VECTOR(bnet)
%
% This function receives a bnet and returns an array containing all
% its parameters following the same ordering as used by the
% constrained learning procedure.
function x = bnet_get_params_as_vector(bnet)
    % get sizes of the bnet
    [nn,ns,np,ni] = bnet_sizes(bnet);
    eqclass = bnet.equiv_class(:);
    % create the parameter array with zeros
    x = zeros(1,sum(ni));
    pos = 1;
    for i = 1:nn
        % for each node, copy to the parameter array the CPT of the node
        % pos and npos keep track of the position in the array x of
        % the values corresponding to the CPTs of node i
        npos = pos + ni(i);
        % verify if it is adjustable
        e = eqclass(i);
        a = struct(bnet.CPD{e});
        if isfield(a,'CPT') && adjustable_CPD(bnet.CPD{e}) && ns(e)>1
            % the correct position to make this copy is determined by the
            % sizes of the preceding CPTs
            x(pos:npos-1) = a.CPT(:)';
        end
        pos = npos;
    end
end
