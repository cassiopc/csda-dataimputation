% % Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% Learn a Bayesian network from missing data
%
% datCell is a cell matrix (samples are in columns, variables in
% rows). bdeu is the ESS (prior strength) for the learning. maxiterEM
% and maxiterSEM are the number of iterations for EM (imputation) and
% for Structure-EM (search for the best graph). typeImp changes how
% the data are imputed within the runs of EM (that is, in the
% E-step). The options are: 
%
% 'DA': impute by sampling and run multiple times EM for imputation
% 'RS': impute by sampling but run it only once for each SEM step
% 'FLDA': impute by sampling running multiple times, but do not run
%         SEM (use a full model). This is used to perform a full model
%         estimation (no independences accounted for) with Data Augmentation
% 'FLE': do not run SEM (use a full model) and impute values by expectation
% 'EE': impute with expectation from the BN and run SEM. This is
%       probably the best performing approach.
% 'EM': impute with the mode from the BN and run SEM. This has to
%       be used in place of EE if the data are categorical with
%       multiple categories (if binary, EE is still ok). 
%
% verb tells whether to be verbose in the output and classe is a
% trick to force a given node to be considered as the class in a
% classification problem. In that case, arcs cannot be used that
% enter the class (this is only reasonable in some particular
% contexts).
%
% Returning values are: the best found Bayesian net in bnetMax, the
% value of the loglik of data (this might be a BDeu score instead
% of likelihood, if a prior smoothing is used). datMax contains the
% imputed data set and topsortMax tells the rearranging of the data
% that was used by the method (the method needs to rearrange the
% order of the variables in the data set because the BNT package
% requires variables to be in topological order with respect to the
% DAG of the Bayesian network -- this is annoying and confusing).
function [bnetMax,maxloglik,datMax,topsortMax] = structureEM(datCell,bdeu,maxiterEM,maxiterSEM,typeImp,verb,classe)
    if nargin < 7
        classe=[];
    end
    if nargin < 6
        verb = 1;
    end
    if nargin < 5
        typeImp = 'E';
    end
    if nargin < 4
        maxiterSEM = 10;
    end
    if nargin < 3
        maxiterEM = 10;
    end
    if nargin < 2
        bdeu = 1;
    end
    fullbnet=0;
    doEM = 1;
    if strcmp(typeImp,'DA')
        typeImp = 'S';
        doEM = -1;
    end
    if strcmp(typeImp,'RS')
        typeImp = 'S';
        doEM = 0;
    end
    if strcmp(typeImp,'FLDA')
        typeImp = 'S';
        doEM = -1;
        maxiterSEM=1;
        fullbnet=1;
    end
    if strcmp(typeImp,'FLE')
        typeImp = 'E';
        doEM = 1;
        maxiterSEM=1;
        fullbnet=1;
    end
    if strcmp(typeImp,'EE')
        typeImp = 'E';
        doEM = 0;
        maxiterSEM=1;
        fullbnet=2;
    end
    if strcmp(typeImp,'EM')
        typeImp = 'M';
        doEM = 0;
        maxiterSEM=1;
        fullbnet=2;
    end
    s = size(datCell);
    if s(1) < 15
        maxpar = 4;
    else
        maxpar = 3;
    end
    if numel(classe) > 0
        allowmat=zeros(s(1),s(1));
        if numel(classe)==numel(allowmat)
            allowmat=classe;
        else
            allowmat(classe,:)=-ones(1,s(1));
            allowmat(:,classe)=ones(1,s(1));
        end
    else
        allowmat=[];
    end
    for i=1:s(1)
        cats{i} = sort(unique(cell2mat(datCell(i,:))));
        ncats{i} = cats{i}(end);
        disp(['var ' num2str(i) ' ncats=' num2str(ncats{i}) ' cats=' num2str(cats{i})]);
    end
    % for i=1:s(1)
    %     for j=1:s(2)
    %         if numel(datCell{i,j}) > 0
    %             datCell{i,j} = find(cats{i}==datCell{i,j});
    %         end
    %     end
    % end
    dag=zeros(s(1));
    if fullbnet==1
        maxpar=s(1);
        for i=1:(s(1)-1)
            for j=(i+1):s(1)
                dag(i,j)=1;
            end
        end
    else
        if fullbnet==2
            allowmat=ones(s(1),s(1));
        end
    end
    bnet = mk_bnet(dag,cell2mat(ncats));
    for i=1:s(1)
        bnet.CPD{i} = tabular_CPD(bnet,i);
    end
    datCompleted=imputation(bnet,datCell,typeImp);
    datMax=datCompleted;
    bnetMax=bnet;
    maxloglik=-inf;
    oldscore=-inf;
    topsort=1:s(1);
    topsortMax=topsort;
    %toporder = 0*bnet.dag;
    %for i=1:s(1)
    %    for j=1:s(1)
    %        if i>j
    %            toporder(i,j)=1;
    %        end
    %    end
    %end
    countr=0;
    if maxiterSEM > 0
        for i=1:maxiterSEM
            %[nn,ns,np,ni] = bnet_sizes(bnet);
            %disp(['BIC=' num2str(sum(ni)*log(s(2))/(2*log(2)))]);
            bnetold = bnet;
            disp(bnet.dag);
            disp(topsort);
            sortold = topsort;
            if verb, disp('call hill-climbing SL method'); end;
            [bnet,bb,topsort] = call_structurelearning(0,[],datCompleted,[],allowmat,-1,16000,verb,6,-1,-1,bdeu,maxpar);
            if doEM~=0 && (numel(topsort)==0 || ~isfield(bb,'score') || bb.score <= oldscore)
                if verb, disp('call exact SL method'); end;
                %algotype=0; % B&B
                if s(1) < 19
                    algotype=1; %dynprog
                else
                    algotype=2;
                end
                gap = oldscore*0.999;
                [bnet,bb,topsort] = call_structurelearning(0,[],datCompleted,[],allowmat,-1,16000,verb,algotype,-1,gap,bdeu,maxpar);
            end
            if numel(topsort)==0 || ~isfield(bb,'score') || bb.score < oldscore
                if isfield(bb,'score') && bb.score < oldscore && countr==0
                    countr=1;
                else
                    break;
                end
            else
                countr=0;
                if verb, disp(['old score=' num2str(oldscore) ' newscore=' num2str(bb.score)]); end;
                oldscore = bb.score;
            end
            if sum(sum(bnetold.dag(sortold,sortold)==bnet.dag(topsort,topsort))) == numel(bnet.dag)
                break;
            end
            if doEM==1
                engine = jtree_inf_engine(bnet);       
                [bnet, LL] = learn_params_em(engine, datCell(topsort,:), maxiterEM,1e-3,1);
                datCompleted(topsort,:)=imputation(bnet,datCell(topsort,:),typeImp);
                if LL(numel(LL))>maxloglik
                    bnetMax=bnet;
                    maxloglik=LL(numel(LL));
                    datMax=datCompleted;
                    topsortMax=topsort;
                end
            else
                datCompleted(topsort,:)=imputation(bnet,datCell(topsort,:),typeImp);
                engine = jtree_inf_engine(bnet);
                [bnet, LL] = learn_params_em(engine, datCompleted(topsort,:), 1,1e-3,1);
                if doEM==-1
                    for k=1:maxiterEM
                        prev = cellmat(datCompleted);
                        datCompleted(topsort,:)=imputation(bnet,datCell(topsort,:),typeImp);
                        dc = cellmat(datCompleted);
                        if any(isnan(prev) | isnan(dc))
                            error('Imputation with mode should not leave any missing on the data');
                        end
                        if all(prev==dc)
                            break;
                        end
                        engine = jtree_inf_engine(bnet);       
                        prevLL = LL(end);
                        [bnet, LL] = learn_params_em(engine, datCompleted(topsort,:), 1,1e-3,1);
                        if LL(end) < prevLL
                            break
                        end
                    end
                end
                bnetMax=bnet;
                maxloglik=LL(numel(LL));
                datMax=datCompleted;
                topsortMax=topsort;
            end
        end
    end
end
%datMax(topsortMax,:)=imputation(bnetMax,datCell(topsortMax,:),typeImp);
