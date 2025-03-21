function [Hn] = comp_hermite_rec(n, x)
%comp_hemrite_dir - Computes the physicist's Hermite polynomial of order n

%Description - Uses recursion (in which the function calls itself) until
%exit conditions are satisfied

%Inputs: n (integer) = any integer greater than or equal to 0
%        x (double scalar) = input value

%Outputs:Hn (double scalar) = evaluation of Hn(x)

%Function dependencies: None

%Author: Edward James, Medical Imaging CDT MRes, 2017

%Define exit conditions for when n = 0 or n = 1
if n == 0
    Hn = 1;
elseif n == 1
    Hn = 2 * x;
elseif n > 1
    %Initiate recursion
    Hn = (2 * x .* comp_hermite_rec(n - 1, x)) - (2 * (n - 1) * comp_hermite_rec(n - 2, x));
end

end
