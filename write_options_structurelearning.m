% % Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% Auxiliary function -- not to be directly invoked
% However one can understand well the available parameters by
% checking the code below. Note that some of these parameters are
% not used by some solvers, so passing by a given parameter does
% not mean it will be used (only structure learning solvers that
% are able to deal with them will).
function write_options_structurelearning(name,varargin)
fd = fopen(name,'w');
args = varargin;
nargs = length(args);
if length(args) > 0
        for i=1:2:(nargs-1)
         if ~isempty(args{i+1})
            switch args{i},
              case 'verbose',      fprintf(fd,'verbose=%d\n',args{i+1});
              case 'maxtreewidth', fprintf(fd,'maxtreewidth=%d\n',args{i+1});
              case 'algorithm',     fprintf(fd,'algorithm=%d\n',args{i+1});
              case 'maxiter',      fprintf(fd,'maxiter=%d\n',args{i+1});
              case 'alpha',     fprintf(fd,'alpha=%d\n',args{i+1});
              case 'mingap',     fprintf(fd,'mingap=%f\n',args{i+1});
              case 'abstarget',     fprintf(fd,'abstarget=%f\n',args{i+1});
              case 'maxmem',     fprintf(fd,'maxmem=%d\n',args{i+1});
              case 'allmaxpar',     fprintf(fd,'globalmaxparents=%d\n',args{i+1});
              case 'maxminutes',     fprintf(fd,'maxminutes=%d\n',args{i+1});
              case 'inputcache',     fprintf(fd,'inputcache=%s\n',args{i+1});
              case 'outputcache',     fprintf(fd,'outputcache=%s\n',args{i+1});
              case 'ilpapprox',     fprintf(fd,'ilpapprox=%d\n',args{i+1});
              case 'ignconstr',     fprintf(fd,'ignconstr=%d\n',args{i+1});
              case 'inibest',     fprintf(fd,'inibest=%d\n',args{i+1});
              case 'maxdata',     fprintf(fd,'maxdata=%d\n',args{i+1});
              case 'node',     fprintf(fd,'node=%d\n',args{i+1});
              case 'threads',     fprintf(fd,'threads=%d\n',args{i+1});
              case 'score',     fprintf(fd,'score=%s\n',args{i+1});
              case 'cacheformat',     fprintf(fd,'cacheformat=%s\n',args{i+1});
              case 'outputfile',     fprintf(fd,'outputfile=%s\n',args{i+1});
              case 'inputfile',     fprintf(fd,'inputfile=%s\n',args{i+1});
            end
        end
      end
end
fclose(fd);
end

