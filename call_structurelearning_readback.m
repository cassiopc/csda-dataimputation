% % Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% Auxiliary function -- not to be directly invoked
function [bnet2,extra,topsort] = call_structurelearning_readback(...
            isdbn,outtemp,nn,ns,verbose,bdeu,reordervars)
    if nargin < 7
        reordervars=1
    end
    fp = fopen(outtemp,'r');
    lin = fgets(fp);
    [msg,er] = ferror(fp);

    % discarding possible initial trash until finding the line with the
    % matrix
    nextra = 0;
    extra = {};
    while ~any(strfind(lin,'matrix::print')) && er == 0
        if verbose ~= 0
            fprintf('%s\n',strtrim(lin));
        end
        nextra = nextra + 1;
        extra{nextra} = strsplit(' ',strtrim(lin));
        lin = fgets(fp);
        [msg,er] = ferror(fp);
    end
    if er ~= 0
        error(['Error in file ' outtemp ]);
    end

    if isdbn
        dag = zeros(nn,nn);
        for j = 1:nn
            for k = 1:nn
                dag(j,k) = fscanf(fp,'%d',1);
            end
        end
        intra = zeros(nn/2);
        inter = zeros(nn/2);
        for j = (nn/2+1):nn
            for k = 1:(nn/2)
                inter(k,j-nn/2) = dag(j,k);
                intra(k,j-nn/2) = dag(j,k+nn/2);
            end
        end
        bnet2 = mk_dbn(intra, inter, ns);
        tok = '';
        % parameters in case of DBN are not being read from the external
        % command but instead shall be learnt later on.
    else
        dag = zeros(nn,nn);
        for j = 1:nn
            for k = 1:nn
                dag(k,j) = fscanf(fp,'%d',1);
            end
        end
        if reordervars
            topsort = topological_sort(dag);
        else
            topsort = 1:nn;
        end
        bnet2 = mk_bnet(dag(topsort,topsort), ns(topsort));
        bnet2.order=1:nn;
        tok = fscanf(fp,'%s',1);
        for jpos = 1:nn
            bnet2.CPD{jpos} = tabular_CPD(bnet2, jpos);
        end
        for j = 1:nn
            if strcmp(tok,'table') ~= 1
                error(['File format error for node ' num2str(j) ]);
            end
            node = fscanf(fp,'%d',1);
            if node+1 ~= j
                error(['File format error for node ' num2str(j) ]);
            end
            tok = fscanf(fp,'%s',1);
            tab = [];
            while strcmp(tok,'table') ~= 1 && strcmp(tok,'score') ~= 1
                tab = [ tab str2num(tok) ];
                tok = fscanf(fp,'%s',1);
            end
            jpos = find(topsort==j);
            ta = struct(bnet2.CPD{jpos});
            if length(ta.CPT(:)) ~= length(tab)
                error(['Table for node ' num2str(j) ' is invalid - perhaps too few data?']);
            end
            parents_size = size(ta.CPT);
            nss = parents_size(end);
            %            parents_size = parents_size((end-1):(-1):(1));
            parents_size = parents_size(1:(end-1));
            parents_size_prod = parents_size;
            parents_size_prod(1) = 1;
            for jj = 2:numel(parents_size)
                parents_size_prod(jj) = parents_size_prod(jj-1) * parents_size(jj-1);
            end
            parents_size = parents_size(end:-1:1);
            indexv = [];
            %disp(['parent_size:' num2str(parents_size)]);
            %disp(['parent_size_prod:' num2str(parents_size_prod)]);
            %disp(['jpos:' num2str(jpos) ' j:' num2str(j)]);
            %disp(['tab:' num2str(tab)]);
            %disp(nss);
            for k = 1:nss
                for jj=1:prod(parents_size)
                    parent_inst = ind2subv(parents_size, jj);
                    parent_inst = parent_inst(end:-1:1);
                    %disp(['inst=' num2str(parent_inst(:)')]);
                    index = 1 + ((parent_inst-1) * parents_size_prod') + (k-1) * prod(parents_size);
                    %disp(['index=' num2str(index)]);
                    %disp(['k=' num2str(k)]);
                    %disp(['jj=' num2str(jj)]);
                    indexv = [indexv, index];
                end
            end
            %disp(['indexv:' num2str(indexv)]);
            %disp(['tab(indexv):' num2str(tab(indexv))]);
            if bdeu < 0
                bnet2.CPD{jpos} = tabular_CPD(bnet2, jpos,'prior_type', 'dirichlet', 'dirichlet_type', 'BDeu', 'dirichlet_weight', -bdeu,'CPT', [tab(indexv)]);
            else
                bnet2.CPD{jpos} = tabular_CPD(bnet2, jpos, [tab(indexv)]);
            end
        end
        if ~bnet_check(bnet2)
            warning('bnet read from result of SL is invalid');
        end
    end
    er = 0;
    while strcmp(tok,'score') ~= 1 && er == 0
        tok = fscanf(fp,'%s',1);
        [msg,er] = ferror(fp);
    end
    if er == 0
        ft = 1;
        lin = fgets(fp);
        [msg,er] = ferror(fp);
        while er == 0
            if ft == 1
                ft = 0;
                lin = sprintf('score %s',strtrim(lin));
            end
            if verbose ~= 0 fprintf('%s\n',strtrim(lin));  end
            nextra = nextra + 1;
            
            extra{nextra} = strsplit(' ',strtrim(lin));
            lin = fgets(fp);
            [msg,er] = ferror(fp); 
        end
            if ft == 1
                lin = sprintf('score %s',strtrim(lin));
            end
        nextra = nextra + 1;
        extra{nextra} = strsplit(' ',strtrim(lin));
        if verbose ~= 0 fprintf('%s\n',strtrim(lin)); end
    end
    fclose(fp);

end
