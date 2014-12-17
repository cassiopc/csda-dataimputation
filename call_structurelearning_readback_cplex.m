% % Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% Auxiliary function -- not to be directly invoked
function [bnet,extra,topsort] = call_structurelearning_readback_cplex(...
            outtemp,nn,ns,data,maxpar,verbose,bdeu,algotype,cachename,reordervars)
    if nargin < 9
        reordervars=1
    end
    cachedata = csvread(cachename);
    dimen = size(cachedata);
    npars = sum(cachedata(:,3:dimen(2))');
    nstat = npars*0;
    for i = 1:length(nstat)
        nstat(i) = cachedata(i,3+cachedata(i,1));
    end
    npars = npars - nstat;
    cachedata = cachedata(find(npars <= maxpar),:);
    [XX,II]=sort(cachedata(:,2));
    cachedata = cachedata(II,:);
    [XX,II]=sort(cachedata(:,1));
    cachedata = cachedata(II,:);

    fp = fopen(outtemp,'r');

    nextra = 0;
    extra = {};
    dag = zeros(nn,nn);

    % discarding possible initial trash until finding the line with the relevant data
    lin = fgets(fp);
    [msg,er] = ferror(fp);
    er = 0;
    while er == 0
        if verbose > 3
            fprintf('READ: %s\n',strtrim(lin));
        end
        if strcmp(sscanf(lin,'%s',2),'Total(root+branch&cut)') == 1
           nextra = nextra + 1;
           extra{nextra} = { 'time', sscanf(lin,'%*s %*s %*s %s',1) };
        end
        if ~isempty(strfind(lin,'maxobj')) && strfind(lin,'maxobj') == 1 
          nextra = nextra + 1;
          extra{nextra} = strsplit(' ',strtrim(lin));
          break;
        end
        if ~isempty(strfind(lin,'absmipgap = ')) && strfind(lin,'absmipgap = ') < 3
          nextra = nextra + 1;
          extra{nextra} = { 'absmingap', sscanf(lin,'%*s %*s %s',1) };
          nextra = nextra + 1;
          extra{nextra} = { 'relmingap', sscanf(lin,'%*s %*s %*s %*s %*s %s',1) };
        end
        lin = fgets(fp);
        [msg,er] = ferror(fp);
    end
    if er ~= 0
        error(['Error in file ' outtemp ]);
    end

    toporder = 1:nn;
    elimorder = 1:nn;
    while 42 == 42
      tex = fscanf(fp,'%s',1);
      fprintf('read type [%s]\n',tex);
      if strcmp(tex,'toporder') == 1
         node = fscanf(fp,'%d',1) + 1;
         fprintf('read node [%d]\n',node);
         pos = fscanf(fp,'%d',1) + 1;
         fprintf('read pos [%d]\n',pos);
         toporder(node) = pos;
      else
       if strcmp(tex,'elimorder') == 1
         node = fscanf(fp,'%d',1) + 1;
         fprintf('read node [%d]\n',node);
         pos = fscanf(fp,'%d',1) + 1;
         fprintf('read pos [%d]\n',pos);
         elimorder(node) = pos;
       else
        if strcmp(tex,'parentnum') == 1
         node = fscanf(fp,'%d',1) + 1;
         fprintf('read node [%d]\n',node);
         pn = fscanf(fp,'%d',1) + 1;
         fprintf('read pn [%d]\n',pn);
         pos = find(cachedata(:,1)==node-1);
         dag(:,node) = cachedata(pos(pn),3:(nn+2));
        else
         break;
        end
       end
      end
    end

    if reordervars
            topsort = topological_sort(dag);
    else
            topsort = 1:nn;
    end
    bnet = mk_bnet(dag(topsort,topsort), ns(topsort),'discrete', [1:nn]);
    bnet.order=topsort;
    for i=1:nn
      bnet.CPD{i} = tabular_CPD(bnet, i,'prior_type', 'dirichlet', 'dirichlet_type', 'BDeu', 'dirichlet_weight', bdeu);
%      bnet.CPD{i} = tabular_CPD(bnet, i, 'CPT', 'unif', bdeu);
    end
    bnet = learn_params(bnet, data(topsort,:));
    if ~bnet_check(bnet)
      fprintf('BNET_CHECK FAILED\n');
    end
    fclose(fp);

end

