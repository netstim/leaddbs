function S = ea_add_StimVector_to_S(S, StimVector,side)

    % Convert OSS-DBS stim vectors to Lead-DBS

    % Stim vector has only cathodic currents (e.g. [-1.0,-0.5,-0.5,-1.0])
    StimVector_perc = -100 * StimVector ./ sum(StimVector);
    S.amplitude{side+1}(:) = 0.0;
    S.amplitude{side+1}(1) = -1.0 * sum(StimVector);

    % anodic case grounding always on
    if side == 0
        S.Rs1.case.perc = 100;
        S.Rs1.case.pol = 2;
    else
        S.Ls1.case.perc = 100;
        S.Ls1.case.pol = 2;
    end

    for cnt=1:length(StimVector)

        if side == 0
            S.Rs1.(['k',num2str(cnt)]).perc = StimVector_perc(cnt);
            if StimVector_perc(cnt) ~= 0
                S.Rs1.(['k',num2str(cnt)]).pol = 1;
                S.activecontacts{1, 1}(1) = 1;
            else
                S.Rs1.(['k',num2str(cnt)]).pol = 0;
                S.activecontacts{1, 1}(1) = 0;
            end
        else
            S.Ls1.(['k',num2str(cnt)]).perc = StimVector_perc(cnt);
            if StimVector_perc(cnt) ~= 0
                S.Ls1.(['k',num2str(cnt)]).pol = 1;
                S.activecontacts{1, 2}(1) = 1;
            else
                S.Rs1.(['k',num2str(cnt)]).pol = 0;
                S.activecontacts{1, 2}(1) = 0;
            end
        end

    end

end