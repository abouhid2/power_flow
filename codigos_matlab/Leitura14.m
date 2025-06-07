%% ABERTURA DO ARQUIVO
FID = fopen('ieee14.PWF','r');  % INDICAR QUAL O NOME DO ARQUIVO A SER LIDO(DEVEM ESTAR NA MESMA PASTA)

cont_row = 1;                % conta o número de linhas do arquivo de dados
cont_dbar = 0;               % conta o numero de barras do sistema
cont_dlin = 0;               % conta o número de linhas e transformadores presentes no sistema

%VARIÁVEIS DE BARRA
number_bar = [];             %número da barra
operation = [];              %operação
state = [];                  %estado
type = [];                   %tipo da barra
group_tension = [];          %grupo base de tensão
name_bar = [];               %nome da barra
group_limits = [];           %grupo limite de tensão
tensao = [];                 %tensao
angle_bar = [];              %angulo
pg = [];                     %geração ativa
qg = [];                     %geração reativa
minqg = [];                  %geração reativa mínima
maxqg = [];                  %geração reativa máxima
bar_control = [];            %barra controlada
pc = [];                     %carga ativa
qc = [];                     %carga reativa
shunt = [];                  %componente shunt(capacitor, reator)
area = [];                   %area

DBAR= [];

%VARIÁVEIS DE LINHA
bar_de = [];                 %barra de
abertura_barra_de = [];      %abertura da barra de
operacao = [];               %operação
abertura_barra_para = [];    %abertura da barra para
bar_para = [];               %barra para
circuito = [];               %circuito da barra
estado = [];                 %estado
proprietario = [];           %
resistencia = [];            %resistencia do circuito em %
reatancia = [];              %reatancia do circuito em %
susceptancia = [];           %valor total da susceptancia shunt do circuito em Mvar
tap = [];                    %tap do transformador referido a barra de
mintap = [];                 %valor mínimo para o tap
maxtap = [];                 %valor maximo para o tap
defasagem = [];              %defasagem do transformador defasador
barra_controlada = [];       %o número da barra a ter a tensão controlado pelo trafo com tap
capacidade_normal = [];      %
capacidade_emergencia = [];  %
numero_de_taps = [];         %

DLIN = [];         %


leitura = fgetl(FID);        %leitura da linha    
while(~strcmp(leitura,'FIM'))
    
    if(strcmp(leitura,'DBAR'))
        for ii=1:2
            leitura = fgetl(FID);
            cont_row = cont_row+1;
        end
        while(~strcmp(leitura,'99999'))
            cont_dbar = cont_dbar + 1;        
            DBAR(cont_dbar,1) = str2double(leitura(1,1:5));
            operation(cont_dbar,1) = leitura(1,6);
            state(cont_dbar,1) = leitura(1,7);
            if(~strcmp(leitura(1,8),' '))
                DBAR(cont_dbar,2) = str2double(leitura(1,8));
            else
                DBAR(cont_dbar,1) = 0;
            end
            group_tension(cont_dbar,:) = leitura(1,9:10);
            name_bar(cont_dbar,1:12) = leitura(1,11:22);
            group_limits(cont_dbar,:) = leitura(1,23:24);
            DBAR(cont_dbar,3) = str2double(leitura(1,25:28))/1000;
            DBAR(cont_dbar,4) = str2double(leitura(1,29:32))*pi/180;
            if(~strcmp(leitura(1,33:37),'     '))
                DBAR(cont_dbar,5) = str2double(leitura(1,33:37));
            else
                DBAR(cont_dbar,5) = 0;
            end
            if(~strcmp(leitura(1,38:42),'     '))
                DBAR(cont_dbar,6) = str2double(leitura(1,38:42));
            else
                DBAR(cont_dbar,6) = 0;
            end
            if(~strcmp(leitura(1,43:47),'     '))                       
                DBAR(cont_dbar,7) = str2double(leitura(1,43:47));
            else
                DBAR(cont_dbar,7) = 0;
            end
            if(~strcmp(leitura(1,48:52),'     '))
                DBAR(cont_dbar,8) = str2double(leitura(1,48:52));
            else
                DBAR(cont_dbar,8) = 0;
            end
            bar_control(cont_dbar,:) = leitura(1,53:58);
            if(~strcmp(leitura(1,59:63),'     '))
                DBAR(cont_dbar,9) = str2double(leitura(1,59:63));
            else
                DBAR(cont_dbar,9) = 0;
            end
            if(~strcmp(leitura(1,64:68),'     '))
                DBAR(cont_dbar,10) = str2double(leitura(1,64:68));
            else
                DBAR(cont_dbar,10) = 0;
            end
            if(~strcmp(leitura(1,69:73),'     '))
                
                DBAR(cont_dbar,11) = str2double(leitura(1,69:73));
            else
                DBAR(cont_dbar,11) = 0;
            end
            area(cont_dbar,1) = str2double(leitura(1,74:76));
            
            cont_row = cont_row+1;
            leitura = fgetl(FID);                        
        end
    end
    
    if(strcmp(leitura,'DLIN'))
        for ii=1:2
            leitura = fgetl(FID);
            cont_row = cont_row+1;
        end
        while(~strcmp(leitura,'99999'))
            length(leitura);
            cont_dlin = cont_dlin + 1;        
            
            DLIN(cont_dlin,1) = str2double(leitura(1,1:5));
            
            abertura_barra_de(cont_dlin,1) = leitura(1,6);
            operacao(cont_dlin,1) = leitura(1,8);
            abertura_barra_para(cont_dlin,1) = leitura(1,10);
            
            DLIN(cont_dlin,2) = str2double(leitura(1,11:15));
            
            circuito(cont_dlin,1) = str2double(leitura(1,16:17));
            estado(cont_dlin,1) = leitura(1,18);
            proprietario(cont_dlin,1) = leitura(1,19);
            
            if(~strcmp(leitura(1,21:26),'      '))
                DLIN(cont_dlin,3) = str2double(leitura(1,21:26));
            else
                DLIN(cont_dlin,3) = 0;
            end
            if(~strcmp(leitura(1,27:32),'      '))
                DLIN(cont_dlin,4) = str2double(leitura(1,27:32));
            else
                DLIN(cont_dlin,4) = 0;
            end
            if(length(leitura)>32 && ~strcmp(leitura(1,33:38),'      '))
                DLIN(cont_dlin,5) = str2double(leitura(1,33:38));
            else
                DLIN(cont_dlin,5) = 0;
            end
            if(length(leitura)>38 && ~strcmp(leitura(1,39:43),'     '))
                DLIN(cont_dlin,6) = str2double(leitura(1,39:43));
            else
                DLIN(cont_dlin,6) = 1;
            end
            if(length(leitura)>43 && ~strcmp(leitura(1,44:48),'     '))
                maxtap(cont_dlin,1) = 1/str2double(leitura(1,44:48));
            else
                maxtap(cont_dlin,1) = 0;
            end
            if(length(leitura)>48 && ~strcmp(leitura(1,49:53),'     '))
                mintap(cont_dlin,1) = 1/str2double(leitura(1,49:53));
            else
                mintap(cont_dlin,1) = 0;
            end
            if(length(leitura)>53 && ~strcmp(leitura(1,54:58),'     '))
                DLIN(cont_dlin,7) = str2double(leitura(1,54:58))/180*pi;
            else
                DLIN(cont_dlin,7) = 0;
            end
            if(length(leitura)>58 && ~strcmp(leitura(1,59:64),'      '))
                barra_controlada(cont_dlin,:) = str2double(leitura(1,59:64));
            else
                barra_controlada(cont_dlin,:) = 0;
            end
            if(length(leitura)>64 && ~strcmp(leitura(1,65:68),'    '))
                capacidade_normal(cont_dlin,:) = leitura(1,65:68);
            else
                capacidade_normal(cont_dlin,:) = '    ';
            end
            if(length(leitura)>68 && ~strcmp(leitura(1,69:72),'    '))
                capacidade_emergencia(cont_dlin,:) = leitura(1,69:72);
            else
                capacidade_emergencia(cont_dlin,:) = '    ';
            end            
            if(length(leitura)>72 && ~strcmp(leitura(1,73:74),'  '))
                numero_de_taps(cont_dlin,:) = leitura(1,73:74);
            else
                numero_de_taps(cont_dlin,:) = '  ';
            end           
           
            cont_row = cont_row+1;
            leitura = [];
            leitura = fgetl(FID);                        
        end
    end
    
    leitura = fgetl(FID);
    cont_row = cont_row+1;
end
