% TANGENT_RESIDUAL_CHECK  MATLAB version of the main residual-scaling test.
%
% For a random rank-p projector P and a tangent perturbation
%   E = P*D*(I-P) + (I-P)*D'*P,  scaled so that ||E||_F = ep,
% the identity Ptilde^2 - Ptilde = E^2 holds exactly, so the idempotency
% defect is O(ep^2) and bounded by ep^2/sqrt(2). A generic symmetric
% perturbation of the same norm gives an O(ep) defect.

rng(0);
n = 100; p = 10; trials = 30;
epsilons = logspace(-6, -1, 6);

res_t = zeros(numel(epsilons), 1);
res_n = zeros(numel(epsilons), 1);
I = eye(n);

for k = 1:numel(epsilons)
    ep = epsilons(k);
    rt = zeros(trials, 1); rn = zeros(trials, 1);
    for t = 1:trials
        [Q, ~] = qr(randn(n, p), 0);
        P = Q * Q';

        D = randn(n);
        E = P*D*(I - P) + (I - P)*D'*P;
        E = E * (ep / norm(E, 'fro'));
        Pt = P + E;
        rt(t) = norm(Pt*Pt - Pt, 'fro');

        S = randn(n); S = (S + S') / 2;
        S = S * (ep / norm(S, 'fro'));
        Pn = P + S;
        rn(t) = norm(Pn*Pn - Pn, 'fro');
    end
    res_t(k) = mean(rt);
    res_n(k) = mean(rn);
end

T = table(epsilons', res_t, res_t ./ (epsilons.^2)', res_n, res_n ./ epsilons', ...
    'VariableNames', {'epsilon', 'TangentResidual', 'TangentOverEps2', ...
                      'GenericResidual', 'GenericOverEps'});
format short e
disp(T);
fprintf('Upper bound for TangentOverEps2: 1/sqrt(2) = %.4f\n', 1/sqrt(2));
