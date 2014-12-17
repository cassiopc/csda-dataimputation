% % Copyright 2008 Information Systems Lab, ECSE Department, Rensselaer
% % Polytechnic Institute (RPI). All rights reserved.
% % This work is licensed under a Creative Commons
% % Attribution-Noncommercial-Share Alike 3.0 United States License
% % http://creativecommons.org/licenses/by-nc-sa/3.0/us/
function name = tmpname
    name = strrep(strrep(tempname('.'),'.\',''),'./','');
    name = sprintf('/tmp/%s',name);
end
