coef = (rand(500,1)-0.5)/30;

typePhiList = {'full', 'reduced'};

for parity = [0,1]
    for targetPre = [true, false]
        for iType = 1:length(typePhiList)
            opts.print = 0;
            opts.targetPre = targetPre;
            opts.typePhi = typePhiList{iType};
            opts.method = 'FPI';
            [phi_FPI,out1] = QSP_solver(coef,parity,opts);
            runtime_FPI = out1.time;
            error_FPI = out1.value;
            opts.method = 'INFFT';
            [phi_INFFT,out2] = QSP_solver(coef,parity,opts);
            runtime_INFFT = out2.time;
            error_INFFT = out2.value;
            fprintf('parity: %d, targetPre: %d, typePhi: %s, norm diff: %.4e\n', parity, opts.targetPre, opts.typePhi, norm(phi_FPI-phi_INFFT));
            fprintf('FPI runtime: %.4f, INFFT runtime: %.4f\n', runtime_FPI, runtime_INFFT);
            fprintf('FPI error: %.4e, INFFT error: %.4e\n', error_FPI, error_INFFT);
        end
    end
end
                