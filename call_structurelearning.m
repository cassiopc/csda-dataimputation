% % Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% The arguments (in order) are
% isdbn,dag,data,constr,allowmat,timelim,memlim,verbose,algotype,maxiter,mingap,bdeu,maxpar,cnames,reordervars,twbound
%
% - bnet: original bnet to be learned (which is not used in fact during the
% learning, but it is just to have a copy of the bnet object to help the
% coding of the function). That is, it may be a bnet with a dag filled
% all with zeros.
%
% - data: is a cell matrix (one line per variable, each column is a
% complete sample/slice). Note that they must be indexed from 1 to
% maxcategory (the function write_file_structurelearning will already do
% the decrement because the slep code is expecting from 0 to maxcategory-1).
% if expectations are to be used, then each position with expected values
% must be an array with the probabilities of the respective categories. It
% looks like this (example with 3 variables per line, first and second
% are binary and the third is ternary):
% 0 1 0
% 1 0 2
% [0.1,0.9] 1 [0.5,0.3,0.2]
% 0 1 [0,0.2,0.8]
% [0.5,0.5] 1 2
%
% - constr (list of parameter or structural constraints). I guess you
% are not going to use them at this moment, but if you want, the format
% must be the same as before for the structure learning code. There is a
% matlab code to generate random constraints also inside the attached.
%
% - allowmat (matrix of hard structure constraints). This is that squared
% matrix with a 1 where an arc may NOT be used, and -1 where an arc is
% enforced. Each line/column is related to a variable. A number 1 in the
% position (i,j) means that the variable i cannot have j as child (which is
% the same as saying that j cannot have i as parent). A number -1 means
% that i must have j as child. Note that the behavior is DIFFERENT when
% using with DBNs.
% In this case, allowmat must be a cell array with exactly two positions.
% The first must contain the matrix with -1,0,1 for the intra-matrix of the
% DBN, and the second element of the cell array must contain the matrix
% with the inter-matrix hard arc constraints.
%
% - timelimit (in minutes). After that, the code stop and return the
% best found network so far. A warning is printed in this case.
%
% - memlimit (in megabytes). The current code does not stop because of
% memory limit. It is smarter in the sense that it keeps running without
% using more memory than that, but the branch and bound search is mixed
% with a backtracking search to guarantee that the memory consumption
% does not exceed that value. In the case of the dynamic programming code,
% then the code stops in case the memory limit is not enough to solve the
% problem. A warning is given in this case.
%
% - verbose (1 to be, 0 to be quiet). Note that even with verbose, the
% matlab does not show the output of the slep.exe code until the end (I
% think because it caches it). Still, with this option on you have some
% additional information in the end.
% 
% algotype: 0 is the default, which means the B&B. 1 indicates to use the
% dynamic programming code, and 2 to use the guessordering code. Code 5
% enhances the guessordering code by running a localsearch after every
% guessordering processing to find the best structure within that order.
%
% maxiter: maximum number of iterations that the code will run. It is
% specially important when running guessordering, because it defines how
% many times to guess an order, so it totally influences the running time.
%
% mingap: used by the B&B procedure to stop when the best current solution
% found so far is already at most mingap (in percentage) always from the
% best possible one. Note that the current solution may even be already the
% best, but the error is related to the best current certificate of that
% solution, which guarantees the maximum error (although it may have
% already no error).
%
% bdeu: weight of the prior, that is, the Equivalent Sample Size,
% a.k.a. the sum of the weights of the Dirichlet prior
%
% maxpar: bound on the maximum number of parents each node may have
%
% cnames: cell array with one string name per variable in the data
% set, used to make the output more readable (when for example the
% graph of the network is drawn)
%
% reordervars: if non-zero, tells the method to reoder the
% variables in topological order when constructing the bnet that is
% returned. This is essential if the bnet is to be later used
% within the BNT package, because BNT demands variables to be in
% topological order. the output parameter topsort gives the order
% in which the variables were put after this transformation. For
% example, if one wants just to check the KL divergence of the
% original network and the one generated here, they have to put
% reordervars equal to zero, and later call bnet_kl_diffgraph. Note
% that to use the old dataset with the new bnet, the permutation of
% variables have to be taken into account! The new permutation can
% be seen by looking to cnames{topsort}, ie, using topsort to
% reindex cnames (well, that if you have built a cnames array of names).

% _NOTE__ABOUT__RUNNING__WITH__DBNs_
% The code is running by duplicating the lines of the data (so twice as
% many variables are created) such that the next sample/slice is put 
% besides the current one. So a data matrix as [ A, B, C, D] is becoming:
% A B C
% B C D
% (one slice is being discarded). Then it learns the graph of the second
% half of variables such that they may have links coming from the first
% half (this will be the links between slices), and links from the same
% second half (which are later used as the links inside a frame).

function [bnet2,bbresult,topsort] = call_structurelearning(...
    isdbn,dag,data,constr,allowmat,...
    timelim,memlim,verbose,algotype,maxiter,mingap,...
    bdeu,maxpar,cnames,reordervars,twbound,cachename,slepexepath,amplexepath,threads,wwhich)
    if nargin < 18 || isempty(slepexepath)
        %slepexepath = '"C:\Documents and Settings\Bill\Desktop\sl\slep.exe"';
        slepexepath = '/usr/bin/slep';
    end
    if nargin < 19 || isempty(amplexepath)
        amplexepath = '/usr/bin/ampl';
    end
    if nargin < 20
        threads=1;
    end
    if nargin < 21
        wwhich=-1;
    end            
    bnet2 = {};
    bbresult = {};
    topsort=[];
    nn=0;
    if numel(dag) < 1
        [aa1,aa2]=size(data);
        dag = zeros(aa1,aa1);
        if numel(dag) > 0
            nn = numel(dag(1,:));
        end
    end
    if numel(data) > 0 && nn > 0
        ns = zeros(1,nn);
        [v1,v2] = size(data);
        for i = 1:v1
            for j = 1:v2
                if numel(data{i,j}) > 1
                    ns(i) = numel(data{i,j});
                    break;
                end
                if numel(data{i,j}) > 0
                    ns(i) = max(ns(i),data{i,j});
                end
            end
            if ns(i) < 2
                %disp(data(i,:));
                error(['variable (' num2str(i) ') must have at least 2 states -- enforce it']);
                %ns(i) = 2;
            end
        end
        % [nn,ns,np] = bnet_sizes(bnet);
    end
    if nargin < 16
        twbound=-1;
    end
    if nargin < 17
        cachename='';
        assert(nn > 0);
    else
        if exist(cachename,'file')
            cachedata = csvread(cachename);
            nn = max(cachedata(:,1))+1;
            ns = zeros(nn,1);
            for i = 1:nn
                a = find(cachedata(:,1)==i-1);
                ns(i) = cachedata(a(1),2+i);
            end
        end
    end
    if nargin < 15
        reordervars=1;
    end
    if nargin < 14
        cnames={};
    end
    if nargin < 10
        % by default, maxiter is never reached :)
        maxiter = 1000000000;
    end
    if nargin < 13
        maxpar = -1;
    end
    if nargin < 12
        % by default, use MDL/BIC
        bdeu = 0;
        score='bic';
    end
    if bdeu > 0
        % in the code input, the prior is given as a negative number... :)
        score='bdeu';
    end
    if nargin < 11
        % by default, mingap is zero (run until exact solution is found)
        mingap = 0;
    end
    if nargin < 9
        % by default, use the B&B algorithm
        algotype = 0;
    end
    if nargin < 8
        % by default, do not be verbose
        verbose = 0;
    end
    if nargin < 7 || memlim < 0
        % by default, 2GB at most of memory to use
        memlim = 2000;
    end
    if nargin < 6 || timelim < 0
        % default time limit of 60 minute
        timelim = 10;
    end
    iguess=-1e99;
    if mingap < -1
        iguess = mingap;
        mingap = mingap * 1.000001;
    end
    dlim=10000000;
    igncons=0;

    if algotype == 4 && (isempty(cachename) || ~exist(cachename))
        error('In order to run the ILP solver you must have created the score cache priorly');
    end
    if algotype == 3 && isempty(cachename)
        error('In order to build the score cache you must provide the cachename where it will be written');
    end
    if algotype == 10 && isempty(amplexepath)
        error(['In order to build the ampl code you must provide the ' ...
               'filename (use amplexepath) where it will be written']);
    end

    optname=tmpname;
    if verbose
        fprintf('building options file: %s\n',optname);
    end
    if numel(data) > 0
        intemp = tmpname;
    else
        intemp = [];
    end
    outtemp = tmpname;
    algotype1 = algotype;
    if algotype1 == 10
        algotype1 = 4;
    end
    write_options_structurelearning(optname,...
                                    'inputfile',intemp,...
                                    'outputfile',outtemp,...
                                    'outputcache',cachename,...
                                    'inputcache',cachename,...
                                    'algorithm',algotype1,...
                                    'maxmem',memlim,...
                                    'maxminutes',(timelim*threads),...
                                    'inibest',iguess,...
                                    'mingap',mingap,...  
                                    'verbose',verbose,...
                                    'score',score,...
                                    'alpha',bdeu,...
                                    'ignconstr',igncons,...
                                    'maxdata',dlim,...
                                    'node',wwhich,...
                                    'maxiter',maxiter,...
                                    'allmaxpar',maxpar,...
                                    'threads',threads,...
                                    'cacheformat','csv',...
                                    'maxtreewidth',twbound);

    if isdbn
        allowmat1 = zeros(nn,nn);
        for i = 1:(nn/2)
            %        for j = (nn/2+1):nn
            for j = 1:nn
                allowmat1(i,j) = 1;
            end
        end
        if (nargin < 4) || (numel(allowmat) == 0)
            allowmat = allowmat1;
        else
            if ~iscell(allowmat) || length(allowmat) ~= 2
                error('For DBN, allowmat must be a cell of two elements. The first is a matrix with the constraints for intra-matrix. The second for inter-matrix');
            end
            for i = 1:2
                [naa nab] = size(allowmat{i});
                if naa ~= nn/2 || nab ~= nn/2
                    error(sprintf('Size of alowmat{%d} is wrong',i));
                end
            end
            for i = (nn/2+1):nn
                for j = 1:(nn/2)
                    allowmat1(i,j) = allowmat{2}(j,i-nn/2);
                    allowmat1(i,j+nn/2) = allowmat{1}(j,i-nn/2);
                end
            end
            allowmat = allowmat1;
        end 
        [ n m ] = size(data);
        data = [ data(:,1:(m-1)) ; data(:,2:m) ];
    else
        allowmat = allowmat';
    end

    if numel(data) > 0
        if verbose
            fprintf('writing input file: %s\n',intemp);
        end
        write_file_structurelearning(intemp,dag,constr,allowmat,[],[],data, cnames);
    end
    %disp('pause'); pause;
    errtemp = tmpname;
    err2temp = tmpname;
    cmd = sprintf('%s %s > %s 2> %s',slepexepath,optname,errtemp,err2temp);
    if verbose
        fprintf('calling: %s\n',cmd);
    end
    rcode = system(cmd);
    if algotype == 10
        if rcode == 0
            copyfile(outtemp,amplexepath);
            fprintf('ampl code save in %s\n', amplexepath);
        else
            fprintf('Warning: non-zero code from slep command: %d. Output follows.\n',rcode);
            fileread(errtemp)
            if exist('err2temp') && exist(err2temp)
                fileread(err2temp)
            end
        end
        delete(err2temp);
        delete(errtemp);
        delete(outtemp);
        delete(optname);
        return
    else
        if rcode == 0 && algotype == 4
            cmd = sprintf('%s %s > %s 2> %s', amplexepath,outtemp,errtemp,err2temp);
            if verbose
                fprintf('calling: %s\n',cmd);
            end
            rcode2 = system(cmd);
        end
    end
    %// retcode >= 10: no solution found
    %// retcode = 3 or 13: iteration stop
    %// retcode = 4 or 14: memory stop
    %// retcode = 5 or 15: mingap stop
    %// retcode = 6 or 16: time-limit stop
    errmsg = { '', '', 'Stopped because of iteration limit', ...
               'Stopped because of memory limit', ...
               'Stopped because of desired gap was achieved', ...
               'Stopped because of time limit', ...
               'No improvement in the local search', ...
               'ERROR' };

    % verify if a error code is returned
    if algotype ~= 7
        if rcode > 28  % 18??
            fprintf('Warning: unknown return code from external command: %d. Output follows.\n',rcode);
            fileread(errtemp)
            if exist('err2temp') && exist(err2temp)
                fileread(err2temp)
            end
            fileread(outtemp)
        else
            if rcode >= 10
                if rcode >= 13
                    fprintf('message: %s\n', errmsg{rcode-10});
                end
                fprintf('Warning: no better solution found.\n');
                %fileread(errtemp)
                %fileread(err2temp)
                %fileread(outtemp)
            else
                if rcode >= 3
                    fprintf('%s\n.', errmsg{rcode});
                end
            end
            if algotype < 8 && algotype ~= 3
                if algotype == 4
                    [bnet2,extra,topsort] = call_structurelearning_readback_cplex(...
                        errtemp,nn,ns,data,maxpar,verbose,bdeu,algotype,cachename,reordervars);
                else
                    [bnet2,extra,topsort] = call_structurelearning_readback(...
                        isdbn,outtemp,nn,ns,verbose,bdeu,reordervars);
                end
                bbresult = GetBBResult(extra);
            else
                bbresult = outtemp;
                bnet2 = {};
                topsort = [];
            end
            if exist('intemp') && exist(intemp)
                delete(intemp);
            end
        end
        if rcode < 3
            if algotype < 8
                delete(outtemp);
            end
            if exist('err2temp') && exist(err2temp)
                delete(err2temp);
            end
            delete(errtemp);
            delete(optname);
        end
    end
end


function bb_result = GetBBResult(raw_result)
    for i = 1 : length(raw_result)
        val = str2num(raw_result{i}{2});
        switch(raw_result{i}{1})
          case 'time'
            bb_result.time = val;
          case 'absmingap'
            bb_result.absgap = val;
          case 'relmingap'
            bb_result.gap = val;
          case 'maxobj'
            bb_result.score = val;
          case 'cachesize'
            bb_result.cachesize = val;
          case 'cachemem'
            bb_result.cachemem  = val;
          case 'totalsearchspace'
            bb_result.totalsearchspace = val;
          case 'reducedsearchspace'
            bb_result.reducedsearchspace = val;
          case 'maxmem'
            bb_result.maxmem = val;
          case 'iter'
            bb_result.iter = val;
          case 'bbscore'
            bb_result.bbscore = val;
          case 'score'
            bb_result.score = val;
          case 'gap'
            bb_result.gap = val;
        end    
    end

end

