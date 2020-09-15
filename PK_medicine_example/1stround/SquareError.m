function sq_error = SquareError(Gfilter, x, y)
%the square error between yhat = Gfilter*x and y
%G is Ny*Nx, 
%x is Nx*nn, 
%y is Ny*nn
vec_error = Gfilter*x - y;    
sq_error = mean(sum(vec_error.^2, 1));
end