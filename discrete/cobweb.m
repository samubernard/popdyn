function [xx,yy] = cobweb(x)
% COBWEB cobweb map from the orbit X
%   [xx,yy] = cobweb(x) is the x and y component of a cobweb
%   map for the orbit X. X must be a row vector. The cobweb can be 
%   plotted with 
%                   plot(xx,yy)


xx = repmat(x,2,1);
n = length(x);
xx = reshape(xx,1,2*n);

yy = xx(2:end-1);
xx = xx(1:end-2);