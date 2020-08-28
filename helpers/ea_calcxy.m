function [x, y] = ea_calcxy(head, tail)
% Calculate X and Y norm vector based on the head and tail of the lead
%
% The X and Y axis orientation will be the same as in the null model of the
% directed lead

trajvector = diff([head; tail]);
normtrajvector = trajvector/norm(trajvector);
y = [0 normtrajvector(3) -normtrajvector(2)];
x = -cross(normtrajvector,y);
y = y/norm(y);
x = x/norm(x);
