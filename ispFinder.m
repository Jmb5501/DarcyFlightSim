% get parameters

params = getParams();

isp = zeros(20, 20);
cdfactor = logspace(0.25, 3, 20);
throatfactor = logspace(0.25, 3, 20);
pararams = params;

for i = [1:20]
    pararams.CDA_f = params.CDA_f * cdfactor(i);
    pararams.CDA_o = params.CDA_o * cdfactor(i);
    for j = [1: 20]
        pararams.A_throat = params.A_throat * throatfactor(j);
        try
           [~, ~, ~, ~, isp(i, j)] = runSim(pararams, 0);
        catch
            isp(i, j) = 0;
        end
    end
end


    
