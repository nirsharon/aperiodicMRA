function [val, grad] = costFunctionV2(x, rho, mu, P, lambda, sigma)
%
% cost function that includes gradient for Matlab's solver
% rho and x are the variables
% mu, lambda, sigma, and  P are the parameters


n = length(x); x = x(:);
Z = repmat(x,1,n);
P_rho = 0;
for j=1:n
    Z(:,j) = circshift(Z(:,j),j);
    P_rho = P_rho + rho(j)*circshift(eye(n),j);
end
% function evaluation
term1 = norm(Z*rho-mu(:))^2;   
term2 = norm(triu(Z*diag(rho)*Z' -P + sigma^2*eye(n)),'fro')^2;
val = term1 + lambda*term2;
% gradient calculation
rho_part1 =  2*Z'*(Z*rho-mu(:));   
x_part1   = 2*P_rho'*(P_rho*x-mu(:));  
% part 2
rho_part2 = 0;
x_part2   = 0;
pp = rho(end:-1:1);
xx = x(end:-1:1);
val_mat = Z*diag(rho)*Z'-P +sigma^2*eye(n);
for j=1:n
    for m=j:n
         rho_part2 = rho_part2 + ...
             2*val_mat(j,m)*(circshift(xx,j-1).*circshift(xx,m-1))';
           x_part2 = x_part2 + ...
               2*val_mat(j,m)*(circshift(pp,j-1).*circshift(x,j-m)+circshift(pp,m-1).*circshift(x,m-j));
    end
end
% unify two parts of cost function
rho_part  = rho_part1 + lambda*rho_part2(:);
x_part    = x_part1 + lambda*x_part2(:);
% unify to one vector
grad = [x_part(:) ; rho_part(:)];
end

