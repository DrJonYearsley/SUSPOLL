function [value, isterminal, direction] = myEvent(~, Tth, maxTth)
% Function to detect when thorax temperature exceed maxTth
% and then send a signal t stop the integration

value = abs(maxTth - Tth(1)); % Even occurs when Tth == maxTth
isterminal = 1;               % Stop the integration
direction  = 0;

end