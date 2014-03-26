function vecint_numberBase = number_representation(int_number, int_digit, int_base)
	% 輸入三個正整數: number dimension baase
	% Eg: number_representation(999, 4, 7) 想把999表示成 999 = 2*7^3 + 6*7^2 + 2*7 + 5
	% 並回傳 [5 2 6 2]
	% 但由於 matlab 矩陣的 index 都是從 1 開始算起 所以為了方便使用回傳 [6 3 7 3]
	vecint_numberBase = zeros(1, int_digit);
	for j = 1:int_digit
		Remainder = mod(int_number, int_base); % 餘數
		vecint_numberBase(j) = Remainder + 1;
		int_number = (int_number - Remainder) / int_base;
	end
end