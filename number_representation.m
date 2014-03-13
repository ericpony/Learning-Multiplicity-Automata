function vecint_NumberBase = number_representation(int_Number, int_Power, int_Base)
% ��J�T�ӥ����: number ���� Base
% Eg: number_representation(999, 4, 7) �Q��999��ܦ� 999 = 2*7^3 + 6*7^2 + 2*7 + 5
% �æ^�� [2 6 2 5]
% ���ѩ�matlab�x�}��index���O�q1�}�l��_ �ҥH���F��K�ϥΦ^��[3 7 3 6]
vecint_NumberBase = zeros(1, int_Power);
for j = 1:int_Power
    Remainder = mod(int_Number, int_Base); % �l��
    vecint_NumberBase(j) = Remainder + 1;
    int_Number = (int_Number - Remainder) / int_Base;
end

end