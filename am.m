function am()
global mtxdb_TargetPolynomial % Target polynomial
global int_VariableNum % Target polynomail的變數個數
global int_TermNum % Taget polynomial的項數
global int_SymbolNum % 每個變數之取樣點個數
global int_Rank % 紀錄演算法目前的階數
global vecsymdb_Symbol

clc; % 清除 History Command
tic;

%============每個變數之取樣點===============
vecsymdb_Symbol = sym(-5:5); % 每個變數之取樣點
int_SymbolNum = size(vecsymdb_Symbol, 2); % 一個變數的取樣點個數
if(int_SymbolNum < 1)
  	error('Invalid Symbol Vector!');
end

fprintf('取樣點:');
display(vecsymdb_Symbol);
%==========================================

%=====計算Lagrange interpolation basis=====
global  mtxsymdb_LagrangeBasis
mtxsymdb_LagrangeBasis = sym(zeros(int_SymbolNum, int_SymbolNum));
for i = 1:int_SymbolNum
    vecdbsym_Temp = sym(1);
    for j = 1:int_SymbolNum
        if(i ~= j)
            vecdbsym_Temp = polymul(vecdbsym_Temp, [1 -vecsymdb_Symbol(j)]) / (vecsymdb_Symbol(i) - vecsymdb_Symbol(j));
        end
    end
    mtxsymdb_LagrangeBasis(i, :) = vecdbsym_Temp;
end
%==========================================

%===========Target polynomial==============
load input.txt; % 從檔案載入矩陣
mtxdb_TargetPolynomial = input(:, :);
int_VariableNum = size(mtxdb_TargetPolynomial, 2) - 1; % 變數個數
int_TermNum = size(mtxdb_TargetPolynomial, 1);

if(int_VariableNum < 1 || int_TermNum < 1)
  	error('Invalid target function!');
end
%==========================================

%==================STEP1===================
% 主要目的:Initialization
cellvecsymdb_SymbolX = {[]; choose_nonezero_string()};
cellvecsymdb_SymbolY = cellvecsymdb_SymbolX;
int_Rank = 2; % 由於Polynomial Learning有些特殊，需要從int_Rank=2開始

vecsymdb_AcceptingState = member_query(cellvecsymdb_SymbolX{1});
mtxdb_Hankel = sym([1 0; 0 0]); % 空字串預設的函數值為1

% 初始化arydb_MemberQuery。為了避免重複Member Querry，使用一個arydb_MemberQuery把已經問過的值保存起來。
arydb_MemberQuery = sym(zeros(int_SymbolNum, 1, 1));

%===========Algorithm Begin===========
while int_Rank < 100
    % Initialize
    fprintf('演算法目前跑到 %d 階\n', int_Rank);
    vecsymdb_AcceptingState(int_Rank) = member_query(cellvecsymdb_SymbolX{int_Rank});  
    arysymdb_SymbolWeighting = sym(zeros(int_SymbolNum, int_Rank, int_Rank));
    
    for i = 1:int_Rank-1 % Extend Hankel matrix
	  	mtxdb_Hankel(i, int_Rank) = member_query([cellvecsymdb_SymbolX{i} cellvecsymdb_SymbolY{int_Rank}]);
	    mtxdb_Hankel(int_Rank, i) = member_query([cellvecsymdb_SymbolX{int_Rank} cellvecsymdb_SymbolY{i}]); 
    end
    mtxdb_Hankel(int_Rank, int_Rank) = member_query([cellvecsymdb_SymbolX{int_Rank} cellvecsymdb_SymbolY{int_Rank}]);

%===========STEP2================
% 主要目的:計算每個symbol的weighting matrix (cellmtxdb_SymbolWeighting)
    mtxdb_InverseHankel = inv(mtxdb_Hankel);
    for i = 1:int_SymbolNum
        for j = 1:int_Rank-1
			arydb_MemberQuery(i, int_Rank, j) = member_query([cellvecsymdb_SymbolX{int_Rank} vecsymdb_Symbol(i) cellvecsymdb_SymbolY{j}]);
			arydb_MemberQuery(i, j, int_Rank) = member_query([cellvecsymdb_SymbolX{j} vecsymdb_Symbol(i) cellvecsymdb_SymbolY{int_Rank}]);
        end
		arydb_MemberQuery(i, int_Rank, int_Rank) = member_query([cellvecsymdb_SymbolX{int_Rank} vecsymdb_Symbol(i) cellvecsymdb_SymbolY{int_Rank}]);
        arysymdb_SymbolWeighting(i, :, :) = reshape(reshape(arydb_MemberQuery(i, :, :), int_Rank, int_Rank) * mtxdb_InverseHankel, 1, int_Rank, int_Rank);
    end
%================================

%===========STEP3================
    vecint_IndexCounterExample = equivalent_querry(arysymdb_SymbolWeighting, vecsymdb_AcceptingState);
    if(size(vecint_IndexCounterExample, 2) == 0) 
        break; % There is no counter example.
    end
    % There is a counter example.
    i = 1;
    flag = true;
    while(i < int_VariableNum && flag) 
		i = i + 1;
    	
        matrix_Temp = eye(int_Rank);
        for j = 1:i-1
            matrix_Temp = sym(matrix_Temp * reshape(arysymdb_SymbolWeighting(vecint_IndexCounterExample(j), :, :), int_Rank, int_Rank));
        end
        
        % 找出要擴充到 cellvecsymdb_SymbolX 及 cellvecsymdb_SymbolY 的 string
        j = 0;
        while(j < int_Rank && flag) % run through y
	    	j = j + 1; 
            
            acc = sum(matrix_Temp(1, :) .* arydb_MemberQuery(vecint_IndexCounterExample(i), :, j));
	    	if(acc ~= member_query([vecsymdb_Symbol(vecint_IndexCounterExample(1:i)) cellvecsymdb_SymbolY{j}]))
	        	int_Rank = int_Rank + 1;
	        	cellvecsymdb_SymbolX{int_Rank} = vecsymdb_Symbol(vecint_IndexCounterExample(1:i-1));
	        	cellvecsymdb_SymbolY{int_Rank} = [vecsymdb_Symbol(vecint_IndexCounterExample(i)) cellvecsymdb_SymbolY{j}];
	        	flag = false;
	    	end
        end
    end
end
toc;

%===========Algorithm End===========
mtxdb_PolyGuess = am2poly(arysymdb_SymbolWeighting, vecsymdb_AcceptingState);
disp(mtxdb_PolyGuess);
% %======Draw Diagram======
% [x, y] = meshgrid(-3:.3:3, -3:.3:3);
% z = mq_grid(x, y, mtxdb_PolyGuess);
% surfc(x,y,z)
% shading interp
% 
% pause;
end

function vecdbsym_PolyMulResult = polymul(vecdbsym_Poly1, vecdbsym_Poly2)
% 同一個variable的多項式乘法
% Eg: (X^2 + 1) * X = (X^3 +  X)
% => polymul([1 0 1], [1 1]) 會ouput [1 0 1 0]
int_DegPoly1 = size(vecdbsym_Poly1, 2);
int_DegPoly2 = size(vecdbsym_Poly2, 2);
vecdbsym_PolyMulResult = sym(zeros(1, int_DegPoly1 + int_DegPoly2 - 1));

for i = 1:int_DegPoly1
    vecdbsym_PolyMulResult(i : i + int_DegPoly2 - 1) = vecdbsym_PolyMulResult(i : i + int_DegPoly2 - 1) + vecdbsym_Poly1(i) * vecdbsym_Poly2;
end

end

function vecsymdbNonezeroString = choose_nonezero_string()
% 函式目標:希望在 target polynomial 的定義域上上找一個點，其函數值不是零
global int_SymbolNum int_VariableNum
global vecsymdb_Symbol

%vecsymdbNonezeroString = zeros(1, int_VariableNum);
for i = 0:int_SymbolNum ^ int_VariableNum - 1
    vecsymdbNonezeroString = vecsymdb_Symbol(number_representation(i, int_VariableNum, int_SymbolNum));
    if(member_query(vecsymdbNonezeroString) ~= 0)
 		return; % 找到函數值非零的取樣點
    end
end
error('Invalid Target Function'); % 找不到函數值非零的取樣點

end





