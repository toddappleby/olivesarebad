function params = fitMultiVarParams(input,output,flagger,params)
options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',1e6,'MaxIterations',5e5);

size(input)
size(output)
if flagger == 1 %xshift
[params,~,~] = lsqcurvefit(@multiHZNL,[max(max(output)*3) params(1) params(2) params(2) params(2)],input,output,[],[],options);
elseif flagger == 2
[params,~,~] = lsqcurvefit(@multiVarNL,[max(max(output)*3) params(1) params(2) params(1) params(2) params(1) params(2)],input,output,[],[],options);
else
[params,~,~] = lsqcurvefit(@multiGainNL,[max(max(output)*3) params(2) params(1) params(1) params(1)],input,output,[],[],options);
end


end