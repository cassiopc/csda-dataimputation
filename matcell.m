% % Copyright 2014 C. P. de Campos (cassiopc@gmail.com). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
% Auxiliary function -- not to be directly invoked
function cel = matcell(mat)
    s=size(mat);
    cel = cell(s);
    for i=1:s(1)
        for j=1:s(2)
            if ~isnan(mat(i,j))
                cel{i,j}=mat(i,j);
            end
        end
    end
end
