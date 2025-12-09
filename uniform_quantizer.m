function y_q = uniform_quantizer(y, Delta, L)
%UNIFORM_QUANTIZER Uniform mid-rise quantizer with L levels (L odd)
%
%   y_q = uniform_quantizer(y, Delta, L)
%
%   y     : scalar input
%   Delta : quantization step size
%   L     : number of quantization levels (odd)
%
%   输出 y_q 是被量化后的标量

    max_idx = (L - 1)/2;

    % normalized value
    z = y / Delta;

    % clip to dynamic range
    z = max(min(z, max_idx), -max_idx);

    % round to nearest integer
    zq = round(z);

    % reconstruct
    y_q = zq * Delta;
end
