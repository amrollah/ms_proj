%%This function implements Mean Absolute Error to be implemented on the
%%predicted and measured power values. The passed predicted and measured
%%values are expected to have matching timestamps. Function will display
%%error if the lengths of the passed values don't match. For future
%%versions, one can consider supersampling the shorter variable array if
%%the lengths of the variables don't match.

%%function implements the equation : MAE =  E[|(predicted-measured)/measured |]
function res = MAE(predicted, measured)

if (length(predicted) ~= length(measured))
    display('lengths of the passed variables are not equal. Breaking.');
    res = nan;
else
    res = mean(abs(predicted-measured)./abs(measured));
end