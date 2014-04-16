function am()
	global mtxdb_TargetPolynomial 	% Target polynomial
	global int_VariableNum 			% Target polynomail 的變數個數
	global int_TermNum 				% Taget polynomial 的項數
	global int_SampleNum 			% 每個變數之取樣點個數
	global int_Rank 				% 紀錄演算法目前的階數
	global vecsymdb_Sample
	global member_query_cache member_query_cache_index_base
	global interactive
	
	interactive = false
	clc; % 清除 History Command	

	%============每個變數之取樣點===============
	%vecsymdb_Sample = sym(-5:5); % 每個變數之取樣點
	vecsymdb_Sample = sym(-1:1); % 每個變數之取樣點
	int_SampleNum = size(vecsymdb_Sample, 2); % 一個變數的取樣點個數
	if(int_SampleNum < 1)
		error('Invalid Symbol Vector!');
	end

	fprintf('Sample points:');
	display(vecsymdb_Sample);
	%==========================================

	%=====計算Lagrange interpolation basis=====
	fprintf('Computing Lagrange basis...');
	t0 = clock;
	global  mtxsymdb_LagrangeBasis
	mtxsymdb_LagrangeBasis = sym(zeros(int_SampleNum, int_SampleNum));
	for i = 1:int_SampleNum
		vecdbsym_Temp = sym(1);
		for j = 1:int_SampleNum
			if(i ~= j)
				vecdbsym_Temp = polymul(vecdbsym_Temp, [1 -vecsymdb_Sample(j)]) / (vecsymdb_Sample(i) - vecsymdb_Sample(j));
			end
		end
		mtxsymdb_LagrangeBasis(i, :) = vecdbsym_Temp;
	end 
	fprintf('elapsed time %.2f sec.\n', etime(clock, t0));
	%==========================================
	%display(mtxsymdb_LagrangeBasis);
	%pause;
	%===========Target polynomial==============
	%Eg. mtxdb_TargetPolynomial = [ 4 1 2 
	%                              -2 3 4 ]
	% represents polynomial 4.x^1.y^2 - 2.x^3.y^4
	
	if(interactive)
		int_TermNum = input('Number of terms in the target polynomial: ');
		int_VariableNum = input('Number of variables in the target polynomial: ');
	else
		load test.txt; % 從檔案載入矩陣
		mtxdb_TargetPolynomial = test(:, :);		
		int_VariableNum = size(mtxdb_TargetPolynomial, 2) - 1;	% 變數個數
		int_TermNum = size(mtxdb_TargetPolynomial, 1);			% 項數
	end

	if(int_VariableNum < 1 || int_TermNum < 1)
		error('Invalid target function!');
	end
	
	member_query_cache = cell(size(vecsymdb_Sample,2),size(vecsymdb_Sample,2))
	member_query_cache(:)={'x'}
	member_query_cache_index_base = 1
	if(min(double(vecsymdb_Sample))<0)
		member_query_cache_index_base = member_query_cache_index_base - min(double(vecsymdb_Sample))
	end
	
	%==========================================

	%==================STEP1===================
	% 主要目的:Initialization
	cellvecsymdb_SymbolX = {[]; choose_nonezero_string()};
	cellvecsymdb_SymbolY = cellvecsymdb_SymbolX;
	int_Rank = 2; % Polynomial Learning 需要從 int_Rank=2 開始

	vecsymdb_AcceptingState = member_query(cellvecsymdb_SymbolX{1});
	mtxdb_Hankel = sym([1 0; 0 0]); % 空字串預設的函數值為1

	% 為了避免重複做 member query，使用 arydb_MemberQuery 把已經回答過的問題保存起來。
	arydb_MemberQuery = sym(zeros(int_SampleNum, 1, 1));

	%===========Algorithm Begin===========
	fprintf('\nEntering main loop:\n');
	t1 = clock;
	while int_Rank < 100
		% Initialize
		fprintf('Running iteration %02d...', int_Rank);
		t0 = clock;
		vecsymdb_AcceptingState(int_Rank) = member_query(cellvecsymdb_SymbolX{int_Rank});  
		arysymdb_SymbolWeights = sym(zeros(int_SampleNum, int_Rank, int_Rank));
		
		for i = 1:int_Rank-1 % Extend Hankel matrix to (int_Rank x int_Rank) square matrix 
			mtxdb_Hankel(i, int_Rank) = member_query([cellvecsymdb_SymbolX{i} cellvecsymdb_SymbolY{int_Rank}]);
			mtxdb_Hankel(int_Rank, i) = member_query([cellvecsymdb_SymbolX{int_Rank} cellvecsymdb_SymbolY{i}]); 
		end
		mtxdb_Hankel(int_Rank, int_Rank) = member_query([cellvecsymdb_SymbolX{int_Rank} cellvecsymdb_SymbolY{int_Rank}]);

		%===========STEP2================
		% 主要目的:計算每個 symbol 的 weighting matrix (cellmtxdb_SymbolWeighting)
		%mtxdb_InverseHankel = inv(mtxdb_Hankel);
		for i = 1:int_SampleNum
			for j = 1:int_Rank-1
				% arydb_MemberQuery(alphabet[i],x,y) = member_query([x alphabet[i] y])
				arydb_MemberQuery(i, int_Rank, j) = member_query([cellvecsymdb_SymbolX{int_Rank} vecsymdb_Sample(i) cellvecsymdb_SymbolY{j}]);
				arydb_MemberQuery(i, j, int_Rank) = member_query([cellvecsymdb_SymbolX{j} vecsymdb_Sample(i) cellvecsymdb_SymbolY{int_Rank}]);
			end
			arydb_MemberQuery(i, int_Rank, int_Rank) = member_query([cellvecsymdb_SymbolX{int_Rank} vecsymdb_Sample(i) cellvecsymdb_SymbolY{int_Rank}]);
			arysymdb_SymbolWeights(i, :, :) = reshape(reshape(arydb_MemberQuery(i, :, :), int_Rank, int_Rank)/mtxdb_Hankel, 1, int_Rank, int_Rank);
		end
		%================================

		%display(mtxdb_InverseHankel);pause;

		%===========STEP3================
		vecint_IndexCounterExample = equivalent_querry(arysymdb_SymbolWeights, vecsymdb_AcceptingState);
		if(size(vecint_IndexCounterExample, 2) == 0) 
			fprintf('elapsed time %.2f sec.\n', etime(clock, t0));
			break; % There is no counter example. We have found the target polynomial.
		end
		
		% Otherwise, there is a counter example and we have to extend the Hankel matrix.
		i = 1;
		flag = true;
		while(i < int_VariableNum && flag) 
			i = i + 1;
			
			matrix_Temp = eye(int_Rank);
			for j = 1:i-1
				matrix_Temp = sym(matrix_Temp * reshape(arysymdb_SymbolWeights(vecint_IndexCounterExample(j), :, :), int_Rank, int_Rank));
			end
			
			% 找出要加到 cellvecsymdb_SymbolX 及 cellvecsymdb_SymbolY 的 string
			
			j = 0;	
			while(j < int_Rank && flag) % run through y

				j = j + 1; % the jth column of the Hankel matrix
				
				acc = sum(matrix_Temp(1, :) .* arydb_MemberQuery(vecint_IndexCounterExample(i), :, j));
				
				if(acc ~= member_query([vecsymdb_Sample(vecint_IndexCounterExample(1:i)) cellvecsymdb_SymbolY{j}]))
					int_Rank = int_Rank + 1;
					cellvecsymdb_SymbolX{int_Rank} = vecsymdb_Sample(vecint_IndexCounterExample(1:i-1));	% w 
					cellvecsymdb_SymbolY{int_Rank} = [vecsymdb_Sample(vecint_IndexCounterExample(i)) cellvecsymdb_SymbolY{j}];	% \sigma + y
					flag = false;	% break
				end
			end
		end
		fprintf('elapsed time %.2f sec.\n', etime(clock, t0));
	end
	
	fprintf('Main loop takes time %.2f sec.\n', etime(clock, t1));
	
	%===========Algorithm End===========
	fprintf('Running interpolation...');
	t0 = clock;
	mtxdb_PolyGuess = am2poly(arysymdb_SymbolWeights, vecsymdb_AcceptingState);
	display(mtxdb_PolyGuess);
	fprintf('Time for doing interpolation: %.2f sec.\n', etime(clock, t0));
	
	% %======Draw Diagram======
%	[x, y] = meshgrid(-3:.3:3, -3:.3:3);
%	z = mq_grid(x, y, mtxdb_PolyGuess);
%	surfc(x,y,z)
%	shading interp
%	pause;
end

function vecdbsym_PolyMulResult = polymul(vecdbsym_Poly1, vecdbsym_Poly2)
	% 同一個variable的多項式乘法
	% Eg: (X^2 + 1) * X = (X^3 +  X)
	% => polymul([1 0 1], [1 1]) 會輸出 [1 0 1 0]
	int_DegPoly1 = size(vecdbsym_Poly1, 2);
	int_DegPoly2 = size(vecdbsym_Poly2, 2);
	vecdbsym_PolyMulResult = sym(zeros(1, int_DegPoly1 + int_DegPoly2 - 1));

	for i = 1:int_DegPoly1
		vecdbsym_PolyMulResult(i : i + int_DegPoly2 - 1) = vecdbsym_PolyMulResult(i : i + int_DegPoly2 - 1) + vecdbsym_Poly1(i) * vecdbsym_Poly2;
	end
end


function vecsymdb_NonzeroString = choose_nonezero_string()
	% 在 target polynomial 的定義域上找一個點，其函數值不是零
	global int_SampleNum int_VariableNum
	global vecsymdb_Sample

	%vecsymdbNonezeroString = zeros(1, int_VariableNum);
	for i = 0:int_SampleNum ^ int_VariableNum - 1
		vecsymdb_NonzeroString = vecsymdb_Sample(number_representation(i, int_VariableNum, int_SampleNum));
		if(member_query(vecsymdb_NonzeroString) ~= 0)
			return; % 找到函數值非零的取樣點
		end
	end
	error('Invalid Target Function'); % 找不到函數值非零的取樣點
end