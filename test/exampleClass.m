classdef exampleClass
    properties
        H;
    end
    methods
        function obj = exampleClass(H)
            obj.H = H;
        end
        function out = fun(obj, x)
            out = x'* obj.H * x;
        end
        function [c, ceq] = constraint(obj, x)
            c = [];
            ceq = x' * x - 1;
        end
    end
    
end
