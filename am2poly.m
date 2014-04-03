function mtxdb_FinalHypothesis = am2poly(arydb_SymbolWeighting, vecdb_AcceptingState)
	% Input: Automata multiplicity 演算法產生的 weighting 矩陣，以及 accepting States 向量
	% Output: 多變數的 polynomial
	% Goal: 利用 lagrange interpolation 把多項式的矩陣表示方式轉成多項式
	global int_VariableNum 	% Target polynomial 的 variable 個數
	global int_SampleNum 	% 對每個 variable，sample points 的個數
	global int_Rank 		% 演算法目前的階數
	cellvecsymdb_SingleSymbolHypothesis = cell(1, 1); % 紀錄 lagrange interpolation 的結果

	t0 = clock;
	% lagrange interpolation
	for i = 1:int_Rank
		for j = 1:int_Rank
			% Each cellvecsymdb_SingleSymbolHypothesis{j, i} is a univariate polynomial
			cellvecsymdb_SingleSymbolHypothesis{i, j} = lagrange(reshape(arydb_SymbolWeighting(:, i, j), 1, int_SampleNum));
		end
	end
	fprintf('Computing H: %.2f sec.\n', etime(clock, t0));
	
	cellvecsymdb_FinalHypothesis = cell(1, 1);	
	cellvecsymdb_TempSymbolHypothesis = cellvecsymdb_SingleSymbolHypothesis;
	
	t0 = clock;
	% Compute [H(z_1) ... H(z_n)]_1
	for m = 1:int_VariableNum-1
		for i = 1:int_Rank	% ith column
			temp = [];	
			for j = 1:int_Rank	% jth row
				% Add row vector [H(z_1)...H(z_m)]_1 * H(z_{m+1}) to temp
				poly1 = cellvecsymdb_TempSymbolHypothesis{1, j};
				poly2 = cellvecsymdb_SingleSymbolHypothesis{j, i};
				%if (isa(poly1,'double') && poly1 == 0 && isa(poly2,'double') && poly2 == 0)
				%	temp = [temp; 0];
				%else
					temp = [temp; polymul_v2(poly1, poly2)];
				%end
			end			
			% cellvecsymdb_FinalHypothesis{i} is the ith component of the 1st row vector of H(z_1)...H(z_{m+1})
			cellvecsymdb_FinalHypothesis{i} = polyadd(unique(double(temp(:, 2:end)), 'rows'), temp); 
		end
		% H(z_1)...H(z_{m+1})
		cellvecsymdb_TempSymbolHypothesis = cellvecsymdb_FinalHypothesis;
	end
	fprintf('Computing [H(z_1) ... H(z_n)]_1 * r: %.2f sec.\n', etime(clock, t0));
	
	if size(cellvecsymdb_FinalHypothesis{2}(:, :), 2) == 0
		mtxdb_FinalHypothesis = [];
		return;
	end
	
	% 把 symbolic 的數值轉成 double
	mtxdb_FinalHypothesis = zeros(size(cellvecsymdb_FinalHypothesis{2}));
	mtxdb_FinalHypothesis(:, 1) = cellvecsymdb_FinalHypothesis{2}(:, 1) * vecdb_AcceptingState(2);
	mtxdb_FinalHypothesis(:, 2:end) = cellvecsymdb_FinalHypothesis{2}(:, 2:end);

end

