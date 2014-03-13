function vecint_NumberBase = number_representation(int_Number, int_Power, int_Base)
% 輸入三個正整數: number 維度 Base
% Eg: number_representation(999, 4, 7) 想把999表示成 999 = 2*7^3 + 6*7^2 + 2*7 + 5
% 並回傳 [2 6 2 5]
% 但由於matlab矩陣的index都是從1開始算起 所以為了方便使用回傳[3 7 3 6]
vecint_NumberBase = zeros(1, int_Power);
for j = 1:int_Power
    Remainder = mod(int_Number, int_Base); % 餘數
    vecint_NumberBase(j) = Remainder + 1;
    int_Number = (int_Number - Remainder) / int_Base;
end

end