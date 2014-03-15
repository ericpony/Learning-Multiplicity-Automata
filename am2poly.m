function mtxdb_FianllHypothesis = am2poly(arydb_SymbolWeighting, vecdb_AcceptingState)
% Input: Automata multiplicity演算法產生的weighting矩陣，以及accepting States向量
% Output: 多變數的polynomial
% Goal: 利用lagrange interpolation把多項式的矩陣表示方式轉成多項式
global int_VariableNum % Target polynomial的variable個數  
global int_SymbolNum % 對每個variable，sample points的個數
global int_Rank % 演算法目前的階數
cellvecsymdb_SingleSymbolHypothesis = cell(1, 1); % 紀錄lagrange interpolation的結果

% lagrange interpolation
for i = 1:int_Rank
	for j = 1:int_Rank
        cellvecsymdb_SingleSymbolHypothesis{i, j} = lagrange(reshape(arydb_SymbolWeighting(:, i, j), 1, int_SymbolNum));
	end
end

% 利用 lagrange interpolation 找出target polynomia
cellvecsymdb_FianllHypothesis = cell(1, 1);
cellvecsymdb_TempSymbolHypothesis = cellvecsymdb_SingleSymbolHypothesis;
for m = 1:int_VariableNum-1
    for i = 1:int_Rank
        temp = [];
        for j = 1:int_Rank
            temp = [temp; polymul_v2(cellvecsymdb_TempSymbolHypothesis{1, j}, cellvecsymdb_SingleSymbolHypothesis{j, i})];
        end
        cellvecsymdb_FianllHypothesis{i} = polyadd(unique(double(temp(:, 2:end)), 'rows'), temp); 
    end
    cellvecsymdb_TempSymbolHypothesis = cellvecsymdb_FianllHypothesis;
end

% 把symbolic的數值轉成double，方便以後計算
mtxdb_FianllHypothesis = zeros(size(cellvecsymdb_FianllHypothesis{2}));
if size(cellvecsymdb_FianllHypothesis{2}(:, :), 2) == 0
    mtxdb_FianllHypothesis = [];
    return;
end
mtxdb_FianllHypothesis(:, 1) = cellvecsymdb_FianllHypothesis{2}(:, 1) * vecdb_AcceptingState(2);
mtxdb_FianllHypothesis(:, 2:end) = cellvecsymdb_FianllHypothesis{2}(:, 2:end);
end

