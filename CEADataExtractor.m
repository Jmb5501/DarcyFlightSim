%% open Filehandle

fh = fopen("CEADataIsopropanol_Unprocessed.txt");

i = 1;
gam = [];
line = fgetl(fh);
j = 1;
u = 1;
lines = [];
while ~isnumeric(line)
    if(contains(line, "Cp, KJ/(KG)(K)"))
        line = string(line);
        linestrs = split(line, ' '); 
        l = 1;
        for k = 1:length(linestrs)
            if((str2double(linestrs(k)))>0)
                cp(j, l) = str2double(linestrs(k));
                l = l + 1; 
            end
        end
        j = j + 1;
        lines = [lines, line];
    end
    if(contains(line, "GAMMAs"))
        line = string(line);
        linestrs = split(line, ' '); 
        l = 1;
        for k = 1:length(linestrs)
            if((str2double(linestrs(k)))>0)
                gam(u, l) = str2double(linestrs(k));
                l = l + 1; 
            end
        end
        u = u + 1;
    end
    i = i+1;
    line = fgetl(fh);
end

R = cp-(cp./gam);

fclose all;
