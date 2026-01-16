function itr_bpm = computeITR(M, P, trialLengthSec, breakSec)
% computeITR  Compute ITR in bits/min for SSVEP BCI.
%
%   M              : number of targets (classes)
%   P              : accuracy (0..1)
%   trialLengthSec : trial duration (s)
%   breakSec       : inter-trial interval (s)

    if nargin < 4
        breakSec = 0;
    end

    if M <= 1
        itr_bpm = 0;
        return;
    end

    % avoid log2(0)
    epsVal = 1e-12;
    P = max(min(P, 1 - epsVal), epsVal);

    B = log2(M) + P*log2(P) + (1-P)*log2((1-P)/(M-1));  % bits/trial
    T = trialLengthSec + breakSec;                      % sec/trial
    itr_bpm = (B * 60) / T;                             % bits/min
end
