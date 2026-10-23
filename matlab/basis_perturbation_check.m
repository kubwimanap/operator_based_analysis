% BASIS_PERTURBATION_CHECK  Perturb an orthonormal basis and re-orthonormalise.
%
% This script perturbs the basis Q0 -> Q0 + eps*G, re-orthonormalises with QR
% and forms Ptilde = Qtilde*Qtilde'. Because Ptilde is then an exact projector,
% the residual ||Ptilde^2 - Ptilde||_F reported here is at rounding level
% (about 1e-15) for every eps. It illustrates that retracting to the manifold
% removes the idempotency defect; it does NOT test the eps^2/sqrt(2) bound.
% For that test see tangent_residual_check.m.

% Parameters
n = 100; p = 10; trials = 100;
epsilons = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12];

% Fixed baseline projector P*
G0 = randn(n,p);
[Q0,~] = qr(G0,0);
Pstar = Q0*Q0';

avg_residuals = zeros(length(epsilons),1);
avg_back_errors = zeros(length(epsilons),1);

for e = 1:length(epsilons)
    ep = epsilons(e);
    residuals = zeros(trials,1);
    back_errors = zeros(trials,1);

    for t = 1:trials
        % Quadratic perturbation model: perturb the basis
        Qpert = Q0 + ep*randn(n,p);
        [Qtilde,~] = qr(Qpert,0);
        Ptilde = Qtilde*Qtilde';

        % Residual: ||Ptilde^2 - Ptilde||_F
        residuals(t) = norm(Ptilde*Ptilde - Ptilde,'fro');

        % Backward error: ||Ptilde - P*||_F
        back_errors(t) = norm(Ptilde - Pstar,'fro');
    end

    avg_residuals(e) = mean(residuals);
    avg_back_errors(e) = mean(back_errors);
end

% Ratios normalized by eps^2
ratios_res = avg_residuals ./ (epsilons.^2)';
ratios_back = avg_back_errors ./ (epsilons.^2)';

% Display results in scientific notation
format short e
T = table(epsilons', avg_residuals, ratios_res, avg_back_errors, ratios_back, ...
    'VariableNames', {'epsilon','Residual','ResidualRatio','BackError','BackRatio'});
disp(T);
