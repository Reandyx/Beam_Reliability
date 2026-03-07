function [Pf,CI] = failure_probability_mc(delta_samples,delta_allow)

N = length(delta_samples);

failures = delta_samples > delta_allow;

Pf = sum(failures) / N;

CI = 1.96 * sqrt( Pf * (1 - Pf) / N );

end