function mtxdb_FianllHypothesis = am2poly(arydb_SymbolWeighting, vecdb_AcceptingState)
% Input: Automata multiplicity�t��k���ͪ�weighting�x�}�A�H��accepting States�V�q
% Output: �h�ܼƪ�polynomial
% Goal: �Q��lagrange interpolation��h�������x�}��ܤ覡�ন�h����
global int_VariableNum % Target polynomial��variable�Ӽ�  
global int_SymbolNum % ��C��variable�Asample points���Ӽ�
global int_Rank % �t��k�ثe������
cellvecsymdb_SingleSymbolHypothesis = cell(1, 1); % ����lagrange interpolation�����G

% lagrange interpolation
for i = 1:int_Rank
	for j = 1:int_Rank
        cellvecsymdb_SingleSymbolHypothesis{i, j} = lagrange(reshape(arydb_SymbolWeighting(:, i, j), 1, int_SymbolNum));
	end
end

% �Q�� lagrange interpolation ��Xtarget polynomia
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

% ��symbolic���ƭ��নdouble�A��K�H��p��
mtxdb_FianllHypothesis = zeros(size(cellvecsymdb_FianllHypothesis{2}));
if size(cellvecsymdb_FianllHypothesis{2}(:, :), 2) == 0
    mtxdb_FianllHypothesis = [];
    return;
end
mtxdb_FianllHypothesis(:, 1) = cellvecsymdb_FianllHypothesis{2}(:, 1) * vecdb_AcceptingState(2);
mtxdb_FianllHypothesis(:, 2:end) = cellvecsymdb_FianllHypothesis{2}(:, 2:end);
end

