% % Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% Function imputes the data by using the estimations from the
% Bayesian network. The result is a complete data matrix
% (represented with matlab cells).
%
% Arguments are:
% bnet: the Bayesian net according to BNT package, which will be
% used to perform the inferences.
% dataCell: a matrix of cells with missing values to be
% imputed. The missing values are given by empty cells
% type: one of 'E', 'M' or 'S' meaning to perform imputation using
% respectively the expected value, the mode or a sampled value (all
% computed using the Bayesian net).
% single: if given, must be a vector with the value to be used in
% the imputation of each line (so dimension must be equal to number
% of rows in the dataCell) 
function dataCell=imputation(bnet,dataCell,type,single)
    engine0 = jtree_inf_engine(bnet);
    s=size(dataCell);
    for j=1:s(2)
        u=[];
        for i=1:s(1)
            if numel(dataCell{i,j})==0
                if nargin > 3 && numel(single) == s(1)
                    dataCell{i,j} = single(i);
                else
                    u=[u,i];
                end
            end
        end
        if numel(u) > 0
            evidence=dataCell(:,j);
            if type=='E'
                engine = enter_evidence(engine0, evidence);
                for i=1:length(u)
                    marg=marginal_nodes(engine, u(i));
                    dataCell{u(i),j}=marg.T;
                end
            else
                if type=='M'
                    mpe = calc_mpe(engine0, evidence, 1); 
                    mpe=cell2num(mpe);
                    for i=1:length(u)
                        dataCell{u(i),j}=mpe(u(i));
                    end
                else
                    if type=='S'
                        sample=sample_bnet(bnet,'evidence',evidence);
                        for i=1:length(u)
                            dataCell{u(i),j} = sample{u(i)};
                        end
                    else
                        error('Invalid type for imputation');
                    end
                end
            end    
        end
    end
end
