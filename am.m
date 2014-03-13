function am()
global mtxdb_TargetPolynomial % Target polynomial
global int_VariableNum % Target polynomail���ܼƭӼ�
global int_TermNum % Taget polynomial������
global int_SymbolNum % �C���ܼƤ������I�Ӽ�
global int_Rank % �����t��k�ثe������
global vecsymdb_Symbol

clc; % �M�� History Command
tic;

%============�C���ܼƤ������I===============
vecsymdb_Symbol = sym(-5:5); % �C���ܼƤ������I
int_SymbolNum = size(vecsymdb_Symbol, 2); % �@���ܼƪ������I�Ӽ�
if(int_SymbolNum < 1)
  	error('Invalid Symbol Vector!');
end

fprintf('�����I:');
display(vecsymdb_Symbol);
%==========================================

%=====�p��Lagrange interpolation basis=====
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
load input.txt; % �q�ɮ׸��J�x�}
mtxdb_TargetPolynomial = input(:, :);
int_VariableNum = size(mtxdb_TargetPolynomial, 2) - 1; % �ܼƭӼ�
int_TermNum = size(mtxdb_TargetPolynomial, 1);

if(int_VariableNum < 1 || int_TermNum < 1)
  	error('Invalid target function!');
end
%==========================================

%==================STEP1===================
% �D�n�ت�:Initialization
cellvecsymdb_SymbolX = {[]; choose_nonezero_string()};
cellvecsymdb_SymbolY = cellvecsymdb_SymbolX;
int_Rank = 2; % �ѩ�Polynomial Learning���ǯS��A�ݭn�qint_Rank=2�}�l

vecsymdb_AcceptingState = member_query(cellvecsymdb_SymbolX{1});
mtxdb_Hankel = sym([1 0; 0 0]); % �Ŧr��w�]����ƭȬ�1

% ��l��arydb_MemberQuery�C���F�קK����Member Querry�A�ϥΤ@��arydb_MemberQuery��w�g�ݹL���ȫO�s�_�ӡC
arydb_MemberQuery = sym(zeros(int_SymbolNum, 1, 1));

%===========Algorithm Begin===========
while int_Rank < 100
    % Initialize
    fprintf('�t��k�ثe�]�� %d ��\n', int_Rank);
    vecsymdb_AcceptingState(int_Rank) = member_query(cellvecsymdb_SymbolX{int_Rank});  
    arysymdb_SymbolWeighting = sym(zeros(int_SymbolNum, int_Rank, int_Rank));
    
    for i = 1:int_Rank-1 % Extend Hankel matrix
	  	mtxdb_Hankel(i, int_Rank) = member_query([cellvecsymdb_SymbolX{i} cellvecsymdb_SymbolY{int_Rank}]);
	    mtxdb_Hankel(int_Rank, i) = member_query([cellvecsymdb_SymbolX{int_Rank} cellvecsymdb_SymbolY{i}]); 
    end
    mtxdb_Hankel(int_Rank, int_Rank) = member_query([cellvecsymdb_SymbolX{int_Rank} cellvecsymdb_SymbolY{int_Rank}]);

%===========STEP2================
% �D�n�ت�:�p��C��symbol��weighting matrix (cellmtxdb_SymbolWeighting)
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
        
        % ��X�n�X�R�� cellvecsymdb_SymbolX �� cellvecsymdb_SymbolY �� string
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
% �P�@��variable���h�������k
% Eg: (X^2 + 1) * X = (X^3 +  X)
% => polymul([1 0 1], [1 1]) �|ouput [1 0 1 0]
int_DegPoly1 = size(vecdbsym_Poly1, 2);
int_DegPoly2 = size(vecdbsym_Poly2, 2);
vecdbsym_PolyMulResult = sym(zeros(1, int_DegPoly1 + int_DegPoly2 - 1));

for i = 1:int_DegPoly1
    vecdbsym_PolyMulResult(i : i + int_DegPoly2 - 1) = vecdbsym_PolyMulResult(i : i + int_DegPoly2 - 1) + vecdbsym_Poly1(i) * vecdbsym_Poly2;
end

end

function vecsymdbNonezeroString = choose_nonezero_string()
% �禡�ؼ�:�Ʊ�b target polynomial ���w�q��W�W��@���I�A���ƭȤ��O�s
global int_SymbolNum int_VariableNum
global vecsymdb_Symbol

%vecsymdbNonezeroString = zeros(1, int_VariableNum);
for i = 0:int_SymbolNum ^ int_VariableNum - 1
    vecsymdbNonezeroString = vecsymdb_Symbol(number_representation(i, int_VariableNum, int_SymbolNum));
    if(member_query(vecsymdbNonezeroString) ~= 0)
 		return; % ����ƭȫD�s�������I
    end
end
error('Invalid Target Function'); % �䤣���ƭȫD�s�������I

end





