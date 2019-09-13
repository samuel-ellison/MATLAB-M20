function v = Laplacian_2D(u, h, stencil)
% This function takes the parameters:
    % u : the input matrix
    % h : the grid spacing
    % stencil : the approximation method (5 or 9)
% and approxiates the Laplacian of u. The output is the approxiated
% Laplacian of u: a matrix v. The Laplacian is approximated via either a
% 5-point or 9-point stencil depending on input. 

% get the size of the input matrix u
[rows,cols] = size(u);

% error checks
if rows < 2 || cols < 2
    error('input matrix must be 2 dimensional!');
end
if stencil ~= 5 && stencil ~= 9
    error('stencil identifier must be 5 or 9!');
end
if isreal(h) == 0 || h <= 0 || mod(h, 1) ~= 0
    error('grid spacing must be a positive, real integer!');
end

% initialize the output matrix v
v0 = zeros(rows,cols);

% define some coefficients for the 5 and 9 stencil methods
x5 = 1/(h^2);
x9 = 1/(6*(h^2));

% iterate through matrix dimensions
for r = 1:1:rows
    for c = 1:1:cols
        
        % define the indices of direct neighbors
        N = r - 1;
        S = r + 1;
        E = c + 1;
        W = c - 1;
        
        % adjust indices if on an edge
        if S > rows
            S = 1;
        end
        if N < 1
            N = rows;
        end
        if E > cols
            E = 1;
        end
        if W < 1
            W = cols;
        end
        
        % Fill in the output matrix depending on which stencil
        if stencil == 5
            % 5 point stencil
            v0(r,c) = x5*(u(N,c) + u(S,c) + u(r,W) + u(r,E) - 4*u(r,c));
            
        else
            % 9 point stencil
            v0(r,c) = x9*(u(N,W) + 4*u(N,c) + u(N,E) + 4*u(r,E) + u(S,E) + 4*u(S,c) + u(S,W) + 4*u(r,W) - 20*u(r,c));
            
        end
    end
end

% equate the output matrix to the recently filled matrix
v = v0;

end
