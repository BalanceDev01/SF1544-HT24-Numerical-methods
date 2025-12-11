clc, clear;

% initial värden
k1 = 5400;
k2 = 135000;
k_vec = [k1; k2];
ck = 0.71;
cs = 1;

tol = 1e-6;
max_it = 100;
errors = zeros(max_it, 1);

for iter = 1:max_it
    F = transfer_functions(k_vec(1), k_vec(2), ck, cs);
    DF = Jacobian_transfer_functions(k_vec(1), k_vec(2));
    s = -DF \ F;
    k_vec = k_vec + s;
    errors(iter) = norm(F);

    if norm(s) < tol
        fprintf('Konvergerade efter %d iterationer.\n', iter);
        break;
    end
end

errors = errors(1:iter);

% se härlednings bild
konvergens = zeros(iter-2, 1);
for i = 2:iter-1
    konvergens(i-1) = log(errors(i+1) / errors(i)) / log(errors(i) / errors(i-1));
end

disp('konvergens ordning:');
disp(konvergens);

k1_opt = k_vec(1);
k2_opt = k_vec(2);

fprintf('Optimala k:\n');
fprintf('k1 = %f\n', k1_opt);
fprintf('k2 = %f\n', k2_opt);

% Optimala k:
% k1 = 691.582879
% k2 = 182116.440285

% transfer_functions(k1_opt, k2_opt, ck, cs);
