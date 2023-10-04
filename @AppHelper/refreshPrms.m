function refreshPrms(obj)
% Refresh the Parameters list and their values.

% Copyright 2022 The MathWorks, Inc.

    app = obj.App;

    % 1. Update the edit fields and/or variable tabs to show the latest values
    
    simInVarNames = [app.modelData.simInp.Variables.Name];
    for idx = 1:numel(app.modelData.tunableVariables)
        varName = app.modelData.tunableVariables(idx).QualifiedName;
        topVarName = obj.extractTopVarName(varName);
        topVar = app.modelData.simInp.Variables(simInVarNames == topVarName).Value;
        wgt = findobj(app.PrmsGL, 'Tag', topVarName);
        if ~isempty(wgt)
            wgt.Value = topVar;
        else
            app.Helper.fillVariableEditTab(varName, topVar, 'Update');
        end
    end

    % 2. Make sure the top panels are visible
    obj.showPanels(["model", "params"]);
    
end % refreshPrms
