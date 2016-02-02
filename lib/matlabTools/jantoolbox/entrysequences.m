function y = entrysequences(x, maxentry)
%#mex
% y = entrysequences(x, maxentry) builds sequences of the entries of x
%   x must contain positive integers between 1 and maxentry.
%   If maxentry is omitted, the maximal value of x is used.
%   After the call, y has the same size as x and contains
%   for each i the number of occurances of x(i) up to then
%   (thus, y(i) == length(find(x(1:i-1)==x(i)))+1)
%
%   Example: x = [5 4 9 9 5 5 9 1 9 9 8 9],
%     then   y = [1 1 1 2 2 3 3 1 4 5 1 6]
%
%   Note: Defining maxentry will not at all affect y, but
%   if it is known, it will increase the performance a little.
%   (If maxentry is not known, do not calculate it, this function
%   is c-code and calculates it faster!)
%
% 31.5.2000, Jan Poland.
