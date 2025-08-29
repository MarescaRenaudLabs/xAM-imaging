function var_exists = check_if_exist(WS, expression)
    var_exists = evalin(WS, ['exist(''' expression ''',''var'');']);
end
