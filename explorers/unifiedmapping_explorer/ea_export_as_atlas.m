function ea_export_as_atlas(app, ss_val, ff_val, nm_val)

if ss_val
    sides = fieldnames(app.explorer.spotdrawn);
    for i = 1:length(sides)
        spot{i} = app.explorer.spotdrawn.(sides{i});
    end
else
    spot = {};
end

if ff_val
    disctract = app.explorer.fiberdrawn;
    disctract.info.Connectome = app.explorer.calcsettings.connectome;

    if app.ShowPositiveCheckBox.Value
        disctract.info.PosAmount = sum(app.explorer.stats.pos.shown);
        if numel(app.explorer.stats.pos.shown) == 2
            disctract.info.PosAmountRight = sum(app.explorer.stats.pos.shown(1));
            disctract.info.PosAmountLeft = sum(app.explorer.stats.pos.shown(2));
            disctract.info.PosPercentRight = app.ShowAmountRightEditFieldPositive.Value;
            disctract.info.PosPercentLeft = app.ShowAmountLeftEditFieldPositive.Value;
        end
        disctract.info.PosPercent = app.ShowAmountEditFieldPositive.Value;
    end

    if app.ShowNegativeCheckBox.Value
        disctract.info.NegAmount = sum(app.explorer.stats.neg.shown);
        if numel(app.explorer.stats.neg.shown) == 2
            disctract.info.NegAmountRight = sum(app.explorer.stats.neg.shown(1));
            disctract.info.NegAmountLeft = sum(app.explorer.stats.neg.shown(2));
            disctract.info.NegPercentRight = app.ShowAmountRightEditFieldNegative.Value;
            disctract.info.NegPercentLeft = app.ShowAmountLeftEditFieldNegative.Value;
        end
        disctract.info.NegPercent = app.ShowAmountEditFieldNegative.Value;
    end
else
    disctract = struct();
end

if nm_val
    network = app.draw;
else
    network = {};
end
ea_unifiedmapping_exportatlas(disctract, spot,network, app.explorer.ID, 'Tract', app.explorer.negcolor, app.explorer.poscolor);
end
