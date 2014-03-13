function matrix_Output = polyadd(mtxint_Degree, matrix_Input)
int_DiffTermNum = size(mtxint_Degree, 1);
int_TermNum = size(matrix_Input, 1);
matrix_Output = sym(zeros(1, size(matrix_Input, 2)));

int_Count = 0;
for i = 1:int_DiffTermNum
    Acc = 0;
    for j = 1:int_TermNum
        if(min(mtxint_Degree(i, :) == matrix_Input(j, 2:end)) ~= 0)
            Acc = Acc + matrix_Input(j, 1);
        end
    end
    if(abs(Acc) ~= 0)
        int_Count = int_Count + 1;
        matrix_Output(int_Count, :) = [Acc mtxint_Degree(i, 1:end)];
    end
end

end
