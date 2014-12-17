% % Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% Auxiliary function -- not to be directly invoked
function mat = cellmat(cel,type)
    if nargin < 2, type = 'X'; end;
    s=size(cel);
    mat = zeros(s(1),s(2));
    for i=1:s(1)
        for j=1:s(2)
            nv = numel(cel{i,j});
            if nv > 0
                if nv > 1
                    if type == 'M'
                        [a,k]=max(cel{i,j});
                        mat(i,j)=k;
                    else
                        if type == 'E'
                            mat(i,j) = (1:nv) * cel{i,j};
                        else
                            if type == 'R'
                                mat(i,j) = round((1:nv) * cel{i,j});
                            else
                                mat(i,j) = nan;
                            end
                        end
                    end
                else
                    mat(i,j)=cel{i,j};
                end
            else
                mat(i,j)=nan;
            end
        end
    end
end
