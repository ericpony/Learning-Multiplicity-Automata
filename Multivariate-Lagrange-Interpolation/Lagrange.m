function result = Lagrange()
    load polynomial.txt;
    load samples.txt;    
    result = [samples Lagrange2(polynomial, samples)];
end

function lagrange_polynomial = Lagrange2(polynomial, sample_points)
    % Input: A polynomial with deg = d, #vars = n, #terms = m = (d+n)!/d!n!, and a vector of m sample points
    % Output: A (m+1)xm matrix that gives the Lagrange representation of the polynomial
    
    degree = max(sum(polynomial(:,2:end), 2));
    num_var = size(polynomial, 2) - 1;    
    num_term = prod(degree+1:degree+num_var)/factorial(num_var);
    
    if(num_var ~= size(sample_points, 2))
        error('Each sample needs %d coordinates.', num_var);
    end
    if(num_term ~= size(sample_points, 1))
        error('Input polynomial needs %d sample points.', num_term);
    end    
    
    nominal_matrx = [[2 0 0]' [0 2 0]' [0 0 2]' [1 1 0]' [1 0 1]' [0 1 1]' [1 0 0]' [0 1 0]' [0 0 1]' [0 0 0]']';
    base_matrix = zeros(num_term, num_term);
    for i = 1:num_term
        for j = 1:num_term
            base_matrix(i,j) = poly_eval([1 nominal_matrx(j, :)], sample_points(i, :));
        end
    end    
    denum = det(base_matrix);
    if(denum==0)
        error('Base matrix must be nonsingular!');
    end    
    
    coeff_matrix = zeros(num_term, num_term);
    sample_values = zeros(1, num_term);
    for row = 1:num_term
        A = base_matrix;
        A(row, :) = [];
        sample_values(1, row) = poly_eval(polynomial, sample_points(row, :));
        for col = 1:num_term
            B = A;
            B(:, col) = [];
            coeff_matrix(row, col) = det(B)*(1-2*mod(row + col, 2));
        end
    end
%    lagrange_polynomial = [(sample_values*coeff_matrix/denum)' nominal_matrx];
    lagrange_polynomial = coeff_matrix/denum;
end

function val = poly_eval(polynomial, point)
    val = 0;
    num_term = size(polynomial, 1);
    for k = 1:num_term
         val = val + polynomial(k, 1) * prod(point .^ polynomial(k, 2:end));
    end
end
