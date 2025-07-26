function updateC(C,index,app)
    if index == 1
        Cfield = app.initJacobi;
        Cselectbar = app.Ci_select_line;
        Cselectpt = app.Ci_select_point;
    elseif index == 2
        Cfield = app.targetJacobi;
        Cselectbar = app.Cf_select_line;
        Cselectpt = app.Cf_select_point;
    end

    Cfield.Value = C;
    Cselectbar.XData = [C,C];
    Cselectpt.XData = C;

end