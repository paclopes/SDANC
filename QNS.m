% Anex to the Paper:
% Low Delay and Low Cost Sigma-Delta Adaptive Controller for Active Noise Control
% Paulo Lopes

classdef QNS < handle
    properties
        mem
        yd
        L
        q
        x
    end
    methods
        function obj = QNS(L)
            obj.L = L;
            reset(obj);
        end
        
        function reset(obj)
            global P;
            obj.mem = zeros(P,1);
            obj.yd = 0;
            obj.q = 0;
            obj.x = 0;
        end
        
        function  y = step(obj, x)
            y = stepx(obj, x, true);
        end
               
        function  y = stepx(obj, x, quantization)
            global P M
            % P: order; M: number of levels
            obj.x = x;
                        
            obj.mem(1) = obj.mem(1) + x - obj.yd;
            for i=2:P
                obj.mem(i) = obj.mem(i) + obj.mem(i-1) - obj.yd;
            end
            
            if quantization
                level_i = min(max(round(M/4*(obj.mem(i)/obj.L+2)-1/2),0),M-1);
                y = obj.L*(4/M*(level_i'+1/2)-2);
            else
                y = obj.mem(P);
            end
            obj.q = y - obj.mem(P);
            obj.yd = y;
        end
    end
end