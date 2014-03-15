function symdb_RequestValue = member_query(vecsymdb_RequestPoint)
global mtxdb_TargetPolynomial int_VariableNum int_TermNum 
symdb_RequestValue = sym(0);

if(size(vecsymdb_RequestPoint, 2) == 0)
    symdb_RequestValue = sym(1);
    return;
end

if(size(vecsymdb_RequestPoint, 2) ~= int_VariableNum)
    return;
end

for i = 1:int_TermNum
     symdb_RequestValue = symdb_RequestValue + mtxdb_TargetPolynomial(i, 1) * prod(vecsymdb_RequestPoint .^ mtxdb_TargetPolynomial(i, 2:end));
end
% symdb_RequestValue = sym(sin(10*vecsymdb_RequestPoint(1)));
% if(double(symdb_RequestValue) > 0)
%     symdb_RequestValue = sym(1);
% elseif(double(symdb_RequestValue) < 0)
%     symdb_RequestValue = sym(-1);
% end

end
