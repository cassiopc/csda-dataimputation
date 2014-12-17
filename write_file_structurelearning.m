% % Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% Auxiliary function -- not to be directly invoked
%
% mutualinfo is a matrix N by N with the value of mutual info
% between the given nodes (the matrix shall be symmetric)

% name is the filename on which the data will be written

% dag is the dag of the bnet of interest. Can be left empty: []

% mallow are initial constraints on arcs. Can be left empty: []

% cluster is a matrix 1 by N or N by 1 with the number of the
% cluster of each node. If the info is not available, it is enough
% to set it to zeros(1,N)

% data is the N by M data (N nodes, M instances) as a cell matrix

% constr is the list of constraints, from {1} to {k}, where k is
% the number of constraints. For example, to fix the constraints we
% are working with, it would be:
% constr{1} = 'T(10)' % of course 10 shall be chosen accordingly
% constr{2} = 'G(10,3)' % 10 is the number of nodes that are
%                       % allowed to have more than 3 parents
% constr{3} = 'F(10,2)' % 10 is the number of nodes that are
%                       % allowed to have more than 2 parents+children
% constr{4} = 'M(10,0.03)' % 10 is the number of arcs that are
%                          % allowed to have mutual info less than 0.03
% constr{5} = 'C(10,5)' % 10 is the number of arcs, 5 the cluster
%                       % number of course 10 shall be chosen accordingly

function write_file_structurelearning(name,dag,constr,mallow,mutualinfo,cluster,data,names)
[ n m ] = size(data);
hasdag = 1;
if numel(dag) == 0
  dag = zeros(n,n);
  hasdag = 0;
end
if numel(mallow) == 0
  mallow = zeros(n,n);
end
if numel(mutualinfo) == 0
  mutualinfo = zeros(n,n);
end
hascluster=1;
if numel(cluster) == 0
    hascluster=0;
%  cluster = zeros(1,n);
end
if n ~= numel(dag(1,:))
    error('Data and dag are incompatible');
end
[ an am ] = size(mallow);
if an > 0 && (an ~= n || am ~= n) 
    error('Allow matrix and dag are incompatible');
end
if nargin < 8
    names={};
end

fd = fopen(name,'w');
if numel(names)==0
    fprintf(fd,'%d %d %d 0 %d %d\n\n', m, n, numel(constr), hasdag, hascluster);
    % write the constraints to the file
    if numel(constr) > 0
        for i = 1:numel(constr)
            fprintf(fd,'-999: %s\n', constr{i});
        end
    end
else
    fprintf(fd,'%d %d %d 0 %d %d\n\n', m, n, -1, hasdag, hascluster);
end

% write the matrix of hard arc constraints, completely empty
for i = 1:n
    for j = 1:n
        if an > 0
            fprintf(fd,'%d ', mallow(i,j));
        else
            fprintf(fd,'0 ');
        end
    end
    fprintf(fd,'\n');
end

if hasdag
for i = 1:n
  for j = 1:n
    fprintf(fd,'%d ', dag(j,i));
  end
  fprintf(fd,'\n');
end

fprintf(fd,'\n\n');
if hascluster
    for i = 1:n
        fprintf(fd,'%d ', cluster(i));
    end
    fprintf(fd,'\n\n');
    
    for i = 1:n
        for j = 1:n
            fprintf(fd,'%d ', mutualinfo(j,i));
        end
        fprintf(fd,'\n');
    end
    fprintf(fd,'\n');
end

if numel(names)>0
    for i=1:n
        fprintf(fd,'%s ',strrep(names{i},' ','_'));
    end
    fprintf(fd,'\n\n');
end
for i = 1:m
    for j = 1:n
        kf = numel(data{j,i});
        if kf>1
            if kf<=2
                fprintf(fd,'%1.8f ', data{j,i}(1));
            else
                fprintf(fd,'[');
                for k = 1:kf
                    fprintf(fd,'%1.8f',data{j,i}(k));
                    if k == kf
                        fprintf(fd,'] ');
                    else
                        fprintf(fd,',');
                    end
                end
            end
        else
            if floor(data{j,i}+0.5) ~= data{j,i}
                fprintf(fd,'%1.8f ', data{j,i}-1);
            else
                fprintf(fd,'%d ', data{j,i}-1);
            end
        end
    end
    fprintf(fd,'\n');
end

fclose(fd);

end
