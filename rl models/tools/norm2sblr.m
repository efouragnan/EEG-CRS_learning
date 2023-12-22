function [ out ] = norm2sblr(in)
% transformation from gaussian space to alpha space, as suggested in Daw 2009 Tutorial
% MKW October 2017

out = 2*(logsig(in)-.5); 

end

