% % Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% Function to run multiple imputation methods over a data matrix.
function [ntests, mse, acc, names, impdata, bnets, topsorts, modeorexpected,single ] = testimputation(filename,R,C,cols,marperc,mcarperc,verb,bdeu,classe,tipos)
    if nargin < 10
        tipos=1:100;
    end
    if nargin < 9
        classe=[];
    end
    if nargin < 7
        verb = 1;
    end
    if nargin < 8, bdeu=1; end;
    if numel(filename) > 1
        M = filename;
        clear filename;
    else
        M=csvread(filename,R,C);
    end
    if numel(cols)>0
        M=M(:,cols);
    end
    s=size(M);
    % train matrix has some elements discarded either by MCAR or MAR
    train = 1+M-repmat(min(M),s(1),1);
    % MAR missing data
    if marperc > 0
        for j=1:s(2)
            for i=2:s(1)
                if train(i-1,j)==train(i,j) && rand > (1-marperc)
                    train(i,j)=nan;
                end
            end
        end
    end
    % MCAR missing data
    if mcarperc > 0
        train(rand(1,numel(train)) > (1-mcarperc)) = nan;
    end
    train = train';
    % test matrix is complete
    test = 1+M-repmat(min(M),s(1),1);
    test=test';
    % define test cases from the matrix
    ntests = numel(find(isnan(train) & ~isnan(test)));
    V=test(find(isnan(train) & ~isnan(test)));
    acc = [];
    mse = [];
    impdata = {};
    bnets = {};
    names = {};
    single= {};
    nn = 1;
    % train the BN
    cel=matcell(train);
    if any(tipos<4)
        [bnet,ml,dat,topsort]=structureEM(cel,bdeu,20,20,'E',1,classe);
    end
    % BN followed by picking mode
    if any(tipos==1)
        datT(topsort,:)=imputation(bnet,dat(topsort,:),'M'); 
        R = cellmat(datT,'M');
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];
        impdata{nn} = R;
        bnets{nn}=bnet;
        topsorts{nn}=topsort;
        names{nn} = 'BN_mode';
        modeorexpected{nn} = 'M';
        if verb, disp(names{nn}); end;
        nn = nn + 1;
    end

    % BN followed by picking expected value
    if any(tipos==2)
        R = cellmat(dat,'E');
        impdata{nn} = R;
        modeorexpected{nn} = 'E';
        names{nn} = 'BN_expectedvalue';
        bnets{nn}=bnet;
        topsorts{nn}=topsort;
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];
    end
    
    if any(tipos==3)
        % BN followed by picking round of expected value
        R = cellmat(dat,'R');
        impdata{nn} = R;
        modeorexpected{nn} = 'R';
        names{nn} = 'BN_roundedexpectedvalue';
        bnets{nn}=bnet;
        topsorts{nn}=topsort;
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];
    end
    % pick intermediate value
    % R = repmat((max(train)+min(train))/2,s(1),1);
    % impdata{nn} = R;
    %names{nn} = 'intermediatevalue';
    % if verb, disp(names{nn}); end;
    % nn = nn + 1;
    % % compute MSE
    % VV=R(find(isnan(train) & ~isnan(test)));
    % mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
    % acc = [ acc, sum(V==VV)/ntests ];

    % % pick rounded intermediate value
    % R = repmat(round((max(train)+min(train))/2),s(1),1);
    % impdata{nn} = R;
    %names{nn} = 'roundedintermediatevalue';
    % if verb, disp(names{nn}); end;
    % nn = nn + 1;
    % % compute MSE
    % VV=R(find(isnan(train) & ~isnan(test)));
    % mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
    % acc = [ acc, sum(V==VV)/ntests ];

    % use single mean imputation
    R = train;
    RR = train;
    s=size(R);
    singval=zeros(1,s(1));
    singmode=zeros(1,s(1));
    for i=1:s(1)
        v=0;
        nv=zeros(max(R(i,:)),1);
        for j=1:s(2)
            if ~isnan(R(i,j))
                v = v + R(i,j);
                nv(R(i,j)) = nv(R(i,j)) + 1;
            end
        end
        v = v/sum(nv);
        singval(i)=v;
        [m,mm] = max(nv);
        singmode(i)=mm;
        for j=1:s(2)
            if isnan(R(i,j))
                % mean
                R(i,j) = v;
                % mode
                RR(i,j) = mm;
            end
        end
    end
    if any(tipos==4)
        [bnet,ml,dat,topsort]=structureEM(cel,bdeu,1,1,'EE',1,classe);
        impdata{nn} = R;
        bnets{nn}=bnet;
        single{nn}=singval;
        modeorexpected{nn} = 'E';
        topsorts{nn}=topsort;
        names{nn} = 'singlemean';
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];
    end
    if any(tipos==5)
        % now the rounded version
        [bnet,ml,dat,topsort]=structureEM(cel,bdeu,1,1,'EE',1,classe);
        R = round(R);
        impdata{nn} = R;
        single{nn}=singmode;
        topsorts{nn}=topsort;
        bnets{nn}=bnet;
        modeorexpected{nn} = 'R';
        names{nn} = 'roundedsinglemean';
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];
    end
    if any(tipos==6)
        % now the mode
        [bnet,ml,dat,topsort]=structureEM(cel,bdeu,1,1,'EM',1,classe);
        bnets{nn}=bnet;
        single{nn}=singmode;
        topsorts{nn}=topsort;
        impdata{nn} = RR;
        modeorexpected{nn} = 'M';
        names{nn} = 'singlemode';
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=RR(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];
    end
    if any(tipos==7)
        % train the BN
        cel=matcell(train);
        [bnet,ml,dat,topsort]=structureEM(cel,bdeu,20,20,'M',1,classe);
        % BN followed by picking mode
        R = cellmat(dat,'M');
        impdata{nn} = R;
        modeorexpected{nn} = 'M';
        names{nn} = 'modeBN_mode';
        bnets{nn}=bnet;
        topsorts{nn}=topsort;
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];
    end
    if any(tipos==8)
        % train the BN with Romero & Salmeron idea IMPUTE_BN
        cel=matcell(train);
        [bnet,ml,dat,topsort]=structureEM(cel,bdeu,20,30,'RS',1,classe);
        % BN followed by picking mode
        R = cellmat(dat,'M');
        impdata{nn} = R;
        modeorexpected{nn} = 'M';
        names{nn} = 'romeros';
        bnets{nn}=bnet;
        topsorts{nn}=topsort;
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];    
    end
    if any(tipos==9)
        % train the BN with DA algo (converges to posterior distro if structure was fixed...)
        cel=matcell(train);
        [bnet,ml,dat,topsort]=structureEM(cel,bdeu,50,20,'DA',1,classe);
        % BN followed by picking mode
        R = cellmat(dat,'M');
        impdata{nn} = R;
        modeorexpected{nn} = 'M';
        names{nn} = 'DA';
        bnets{nn}=bnet;
        topsorts{nn}=topsort;
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];    
    end
    if any(tipos==10)
        % train the BN with full (complete graph) DA algo (converges to posterior distro...)
        cel=matcell(train);
        [bnet,ml,dat,topsort]=structureEM(cel,bdeu,20,1,'FLDA',1,classe);
        % BN followed by picking mode
        R = cellmat(dat,'M');
        impdata{nn} = R;
        modeorexpected{nn} = 'M';
        names{nn} = 'FullDA';
        bnets{nn}=bnet;
        topsorts{nn}=topsort;
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];    
    end
    if any(tipos==11)
        % train the BN with full (complete graph) EM algo (converges to MAP...)
        cel=matcell(train);
        [bnet,ml,dat,topsort]=structureEM(cel,bdeu,20,1,'FLE',1,classe);
        % BN followed by picking mode
        R = cellmat(dat,'E');
        impdata{nn} = R;
        modeorexpected{nn} = 'E';
        names{nn} = 'FullEM';
        bnets{nn}=bnet;
        topsorts{nn}=topsort;
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];    
    end
    if any(tipos==12)
        % train the BN with full (complete graph) EM algo followed
        % by mode
        cel=matcell(train);
        [bnet,ml,dat,topsort]=structureEM(cel,bdeu,20,1,'FLM',1,classe);
        % BN followed by picking mode
        datT(topsort,:)=imputation(bnet,dat(topsort,:),'M'); 
        R = cellmat(datT,'M');
        impdata{nn} = R;
        modeorexpected{nn} = 'M';
        names{nn} = 'FullEMMode';
        bnets{nn}=bnet;
        topsorts{nn}=topsort;
        if verb, disp(names{nn}); end;
        nn = nn + 1;
        % compute MSE
        VV=R(find(isnan(train) & ~isnan(test)));
        mse = [ mse, sqrt(sum((V-VV)'*(V-VV))/ntests) ];
        acc = [ acc, sum(V==VV)/ntests ];    
    end
end
