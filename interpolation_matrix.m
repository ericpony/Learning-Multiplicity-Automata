function mtx = interpolation_matrix(cellvecdb_Hypothesis, interpolation)
size_matrix = size(cellvecdb_Hypothesis, 1);
mtx = zeros(size_matrix, size_matrix);
for i = 1:size_matrix
    for j = 1:size_matrix
        poly = cellvecdb_Hypothesis{i,j};
        size_poly = size(poly, 1);
        if(size_poly > 1)
            acc = 0;
            for k = 1:size_poly
                acc = acc + poly(k, 1)*interpolation^poly(k,2);
            end
            mtx(i, j ) = acc;
        else
            mtx(i, j ) = 0;
        end
    end
end

end

