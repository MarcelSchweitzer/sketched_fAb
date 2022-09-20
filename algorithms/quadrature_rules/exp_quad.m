function [z,c] = exp_quad(active_nodes,tol,ell)
% quadrature nodes and weights for exp(z), taken from funm_quad
% contour encloses active_nodes 
% contour truncation tol
% number of nodes ell

% uses midpoint rule on parabolic Hankel contour
aa = max(1,max(real(active_nodes))+5); % SG: changed 1 to 5
bb = 1;
thetas = imag(active_nodes);

ccs = abs((active_nodes - aa - 1i*thetas)./thetas.^2);
cc = min(ccs)/5; % safety parameter
cc = min(cc,0.25);

phi = @(theta) aa + 1i*bb*theta - cc*theta.^2;
thetac = sqrt((aa-log(tol))/cc); % critical theta
theta = linspace(-thetac,thetac,ell);
hh = theta(2)-theta(1);

z = phi(theta); % quad nodes
c = hh/(2i*pi)*exp(z).*(1i*bb - 2*cc*theta); %weights
