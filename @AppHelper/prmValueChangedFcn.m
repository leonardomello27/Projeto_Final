function prmValueChangedFcn(obj, event)
%PRMVALUECHANGEDFCN Callback for updating a tunable variable's value.
%   This method is called when the app user edits the value of a variable
%   listed in the model parameters pane located at the right of the model
%   view. It updates the variable's value in the app.modelData structure.
%   This structure is what the Simulink.SimulationInput uses to retrieve
%   the values of variables that are eventually applied to the model during
%   the simulation.

%   Copyright 2021 The MathWorks, Inc.

    app = obj.App;
    % Parameters are displayed as name-value pairs on the UI.
    % To get the index of a parameter within the grid, we count by 2's.
    tunableVarIdx = 2 * event.Source.Layout.Row;
    varName = app.PrmsGL.Children(tunableVarIdx - 1).Text;
    tvIdx = string({app.modelData.tunableVariables.QualifiedName}) == varName;

    obj.debugTrace([varName, ' : ', ...
        num2str(app.modelData.tunableVariables(tvIdx).Value), ' -> ', ...
        num2str(event.Value)]);

    % update value 
    app.modelData.tunableVariables(tvIdx).Value = event.Value;

    if app.simStatus == slsim.SimulationStatus.Paused
        % keep track of parameter changes when simulation is paused, these
        % will be applied when simulation is resumed.
        modTV = app.modelData.tunableVariables(tvIdx);
        if isempty(app.modelData.modifiedTunableVariables)
            app.modelData.modifiedTunableVariables = modTV;
        else
            app.modelData.modifiedTunableVariables(end+1) = modTV;
        end
    end
end
