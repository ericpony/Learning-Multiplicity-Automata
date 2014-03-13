function vecint_IndexCounterExample = equivalent_querry(arysymdb_SymbolWeighting, vecsymdb_AcceptingState)
global int_SymbolNum int_VariableNum vecsymdb_Symbol int_Rank

vecint_Test = randi(int_SymbolNum ^ int_VariableNum, 100, 1);
for i = 1:100 % 隨機測試100個點，若都滿足..就算對了
	mtxsymdb_Temp = zeros(1, int_Rank);
    mtxsymdb_Temp(1) = 1;
    vecint_IndexCounterExample = number_representation(vecint_Test(i), int_VariableNum, int_SymbolNum);
    for j = 1:int_VariableNum
 		mtxsymdb_Temp = mtxsymdb_Temp * reshape(arysymdb_SymbolWeighting(vecint_IndexCounterExample(j), :, :), int_Rank, int_Rank);
    end
    
    if (mtxsymdb_Temp(1, :) * vecsymdb_AcceptingState' - member_query(vecsymdb_Symbol(vecint_IndexCounterExample))) ~= 0
 		return;
    end
end
vecint_IndexCounterExample = []; % 表示學到的函數在取樣點上皆相等

end
