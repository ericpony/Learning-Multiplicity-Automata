function mtxsymdb_MulResult = polymul_v2(mtxsymdb_Polynomail1, mtxsymdb_Polynomail2)
int_Poly1TermNum = size(mtxsymdb_Polynomail1, 1); % 多項式1的term數
int_Poly2TermNum = size(mtxsymdb_Polynomail2, 1); % 多項式2的term數
if(int_Poly1TermNum == 0 && int_Poly2TermNum == 0)
    error('Invalid Polynomial');
end

int_Poly1VarNum = size(mtxsymdb_Polynomail1, 2) - 1;
int_Poly2VarNum = size(mtxsymdb_Polynomail2, 2) - 1; % 多項式2的變數個數
if(int_Poly1VarNum < 1 || int_Poly2VarNum < 1) % 零多項式
    mtxsymdb_MulResult = []; 
    return;
end

int_Count = 0;
mtxsymdb_MulResult = sym(zeros(int_Poly1TermNum * int_Poly2TermNum, int_Poly1VarNum + 2));
for i = 1:int_Poly1TermNum
    for j = 1:int_Poly2TermNum
        int_Count = int_Count + 1;
        mtxsymdb_MulResult(int_Count, 1) = mtxsymdb_Polynomail1(i, 1)*mtxsymdb_Polynomail2(j, 1);
        mtxsymdb_MulResult(int_Count, 2:end) = [mtxsymdb_Polynomail1(i, 2:end) mtxsymdb_Polynomail2(j, 2:end)];
    end
end

end
