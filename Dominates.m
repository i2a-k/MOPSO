function [Out] = Dominates(x,y)

Out = all(x <= y) && any(x < y);