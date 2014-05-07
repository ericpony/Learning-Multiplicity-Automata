function symdb_RequestValue = member_query(vecsymdb_RequestPoint)


	global mtxdb_TargetPolynomial int_VariableNum int_TermNum interactive member_query_cache member_query_cache_index_base vecsymdb_Sample
	symdb_RequestValue = sym(0);
	if(size(vecsymdb_RequestPoint, 2) ~= int_VariableNum)
		return;
	end
	if(size(vecsymdb_RequestPoint, 2) == 0)
		symdb_RequestValue = sym(1);
		return;
	end
	
	
	indices = num2cell(vecsymdb_RequestPoint + member_query_cache_index_base);
	cached = member_query_cache(indices{:});
	cached = cached{1,1};
	if(cached ~= 'x')
		symdb_RequestValue = cached;
		fprintf(strrep(['Cached value for query [' num2str(double(vecsymdb_RequestPoint), ' %d ') ']'], ']', '] is %d\n'), double(cached));
		return
	end
	
	if(interactive)		
		fprintf('=== Member query for points ===')
		display(vecsymdb_RequestPoint)
		symdb_RequestValue = sym(input('Query result: '))
		member_query_cache(indices{:}) = {symdb_RequestValue};
		return 
	end	
	
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
	fprintf('=== Member query ===')
	display(vecsymdb_RequestPoint)
	fprintf('Return value: %d\n', double(symdb_RequestValue));
	member_query_cache(indices{:}) = {symdb_RequestValue};
end
