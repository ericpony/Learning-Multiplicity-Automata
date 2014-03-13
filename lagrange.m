function mtxsymdb_InterpolationPolynomial = lagrange(vecdb_PointY)
global int_SymbolNum % ��C��variable�Asample points���Ӽ�
global mtxsymdb_LagrangeBasis % lagrange interpolation����
if(sum(abs(vecdb_PointY)) == 0) % �ͪ��s�h����
    mtxsymdb_InterpolationPolynomial = sym(0);
    return;
end

% Begin to interpolate 
mtxsymdb_InterpolationPolynomial = sym(zeros(1, int_SymbolNum));
for i = 1:int_SymbolNum
    mtxsymdb_InterpolationPolynomial = mtxsymdb_InterpolationPolynomial + mtxsymdb_LagrangeBasis(i, :) * vecdb_PointY(i);
end

mtxsymdb_InterpolationPolynomial = [mtxsymdb_InterpolationPolynomial' sym(int_SymbolNum-1:-1:0)'];
mtxsymdb_InterpolationPolynomial(abs(mtxsymdb_InterpolationPolynomial(:,1)) == 0, :) = []; % �R���h�������A�Y�Ƶ���s��Term
end

