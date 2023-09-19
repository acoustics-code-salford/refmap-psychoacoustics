function [y] = LowPass(x, r_s)
% 

k = 3; % Footnote 21 ECMA-418-2:2022
e_i = [0, 1, 1]; % Footnote 21 ECMA-418-2:2022

% Footnote 20 ECMA-418-2:2022
tau = 1/32*6/7;

d = exp(-1/(r_s*tau)); % Section 5.1.4.2 ECMA-418-2:2022

% Feed-backward coefficients, Equation 14 ECMA-418-2:2022
m = 1:k;
a = [1, ((-d).^m).*arrayfun(@(m_) nchoosek(k, m_), m)];

% Feed-forward coefficients, Equation 15 ECMA-418-2:2022
m = 0:k-1;
i = 1:k-1;
b = (((1 - d)^k)./sum(e_i(i + 1).*(d.^i))).*(d.^m).*e_i;

% Recursive filter Equation 13 ECMA-418-2:2022
y = filter(b, a, x);

end