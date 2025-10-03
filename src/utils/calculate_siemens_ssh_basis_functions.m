function [basis, order, degree, dummy] = calculate_siemens_ssh_basis_functions(r, theta, phi, R0, nmax)

N = length(r); % total number of voxels

%% Calculate the total number of basis functions
Nl = 0;
for n = 0:nmax
    Nl = Nl + (n + 1);
end

%% Calculate SSH basis functions (N x Nl x 2)
basis = zeros(N, Nl, 2, 'double'); % [m]
dummy = zeros(Nl, 1, 'logical');
order = zeros(Nl, 1, 'double');
degree = zeros(Nl, 1, 'double');

offset = 0;

for n = 0:nmax
    %----------------------------------------------------------------------
    % Calculate the radius term
    %----------------------------------------------------------------------
    r_term = R0 * (r / R0).^n; % [m]

    %----------------------------------------------------------------------
    % Compute the associated Legendre functions of degree n and order m = 0, 1, ..., n
    %----------------------------------------------------------------------
    p_unnormalized = legendre(n, cos(theta)).'; % N x n+1

    %----------------------------------------------------------------------
    % Calculate the number of basis functions for this order
    %----------------------------------------------------------------------
    nr_coefficients = n + 1;

    for m = 0:n
        %------------------------------------------------------------------
        % Calculate the normalization factor
        %------------------------------------------------------------------
        if m > 0
            normalization_factor = (-1)^m * sqrt((2 * n + 1) * factorial(n - m) / (2 * factorial(n + m)));
        else % m == 0
            normalization_factor = 1;
        end
        %fprintf('(n,m)=(%2d,%2d) norm_factor = %+7.5f\n', n, m, normalization_factor);

        %------------------------------------------------------------------
        % Calculate the cosine and sine basis functions
        %------------------------------------------------------------------
        basis(:, offset + m + 1, 1) = r_term .* cos(m * phi) .* (normalization_factor * p_unnormalized(:,m + 1));
        basis(:, offset + m + 1, 2) = r_term .* sin(m * phi) .* (normalization_factor * p_unnormalized(:,m + 1));

        %------------------------------------------------------------------
        % Record the information of the spherical harmonics
        %------------------------------------------------------------------
        if m == 0
            dummy(offset + m + 1) = true;
        end

        order(offset + m + 1) = n;
        degree(offset + m + 1) = m;
    end

    %----------------------------------------------------------------------
    % Calculate the cumulative number of basis functions
    %----------------------------------------------------------------------
    offset = offset + nr_coefficients;
end

end