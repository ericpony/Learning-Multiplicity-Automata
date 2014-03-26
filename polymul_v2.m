function mtxsymdb_MulResult = polymul_v2(mtxsymdb_Polynomail1, mtxsymdb_Polynomail2)
	% 多項式1:	[ 1/2, 4, 1] 代表 1/2*(x^4)(y) + 3*(x^7)(y^2)
	%			[ 3,   7, 2]
	% 多項式2 必須為 univariate, 且其變數和多項式1的變數不重複
	% 計算結果的 column 數為 mtxsymdb_Polynomail1 的 column 數加一(多了一個變數), row 為多項式的項數
	% 
	int_Poly1TermNum = size(mtxsymdb_Polynomail1, 1); 	% 多項式1的項數
	int_Poly2TermNum = size(mtxsymdb_Polynomail2, 1); 	% 多項式2的項數
	if(int_Poly1TermNum == 0 && int_Poly2TermNum == 0)
		error('Invalid Polynomial');
	end

	int_Poly1VarNum = size(mtxsymdb_Polynomail1, 2) - 1;	% 多項式1的變數個數
	int_Poly2VarNum = size(mtxsymdb_Polynomail2, 2) - 1; 	% 多項式2的變數個數
	
	if(int_Poly1VarNum < 1 || int_Poly2VarNum < 1) 			% 零多項式
		mtxsymdb_MulResult = []; 
		return;
	end

	int_MulResultTermNum = 0;	% 計算結果的項數
	%mtxsymdb_MulResult = sym(zeros(int_Poly1TermNum * int_Poly2TermNum, int_Poly1VarNum + 2));
	mtxsymdb_MulResult = sym(zeros(int_Poly1TermNum * int_Poly2TermNum, int_Poly1VarNum + int_Poly2VarNum + 1));
	
	for i = 1:int_Poly1TermNum
		for j = 1:int_Poly2TermNum
			int_MulResultTermNum = int_MulResultTermNum + 1;
			mtxsymdb_MulResult(int_MulResultTermNum, 1) = mtxsymdb_Polynomail1(i, 1)*mtxsymdb_Polynomail2(j, 1);
			mtxsymdb_MulResult(int_MulResultTermNum, 2:end) = [mtxsymdb_Polynomail1(i, 2:end) mtxsymdb_Polynomail2(j, 2:end)];
		end
	end
end
