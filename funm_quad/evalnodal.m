function [ p ] = evalnodal( x,nodes,subdiag )
% returns prod(subdiag) ./ ( (x-nodes(1))*(x-nodes(2))*... )

if length(nodes) ~= length(subdiag),
    error('EVALNODAL: nodes and subdiag not of the same length!');
end

p = 0*x + 1;
for j = 1:length(nodes),
    p = p * subdiag(j) ./ (x - nodes(j));
end

end