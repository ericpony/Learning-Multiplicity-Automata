function mtxsymdb_InterpolationPolynomial = lagrange(vecdb_PointY)
	global int_SampleNum 			% 對每個 variable，sample points 的個數
	global mtxsymdb_LagrangeBasis 	% lagrange interpolation 的基底
	if(sum(abs(vecdb_PointY)) == 0) % 趨近於零多項式
		mtxsymdb_InterpolationPolynomial = sym(0);
		return;
	end

	% Begin to interpolate 
	mtxsymdb_InterpolationPolynomial = sym(zeros(1, int_SampleNum));
	for i = 1:int_SampleNum
		mtxsymdb_InterpolationPolynomial = mtxsymdb_InterpolationPolynomial + mtxsymdb_LagrangeBasis(i, :) * vecdb_PointY(i);
	end

	mtxsymdb_InterpolationPolynomial = [mtxsymdb_InterpolationPolynomial' sym(int_SampleNum-1:-1:0)'];
	mtxsymdb_InterpolationPolynomial(abs(mtxsymdb_InterpolationPolynomial(:,1)) == 0, :) = []; % 刪掉多項式中，係數等於零的Term
end

