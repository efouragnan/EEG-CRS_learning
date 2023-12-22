function [ alpha_out ] = norm2alpha(alpha)
% transformation from gaussian space to alpha space, as suggested in Daw 2009 Tutorial
% MKW October 2017

alpha_out = logsig(alpha); % same as sigmoid(alpha)

end

