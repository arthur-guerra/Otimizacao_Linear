#Big M - Simplex (Otimização Linear_UFPB)
#Arthur Leandro Guerra Pires

using LinearAlgebra
using PrettyTables

        tableau_inicial =[-3 -5 0; 1 0 4 ; 0 2 12 ; 3 2 18]

        vetor_sinal = ["=", "<=" , "<=" , "<="] 

        println("\nO tableau original do problema é: ")

        #"minimização" ou "maximização"

        tipo_otimização = "maximização"
        
         #__________________________________________________________________________________
        # Dimensões do tableau_inicial

        (n_linhas,n_colunas) = size(tableau_inicial)
        
        #__________________________________________________________________________________
        # Impressões iniciais

        header_inicial = String[]

        for j = 1:(n_colunas)
            push!(header_inicial, "x$j")
        end

        header_inicial[n_colunas, 1] = "b"

        pretty_table(tableau_inicial,header = header_inicial, border_crayon = crayon"yellow",formatters = ft_round(2))

        #__________________________________________________________________________________
        # Criando valor M com os indices da FO

        valor_M = sum(tableau_inicial, dims=2)

        M = valor_M[1,1]

        ### println(M)

        if M <= 0
            M = M *-1000
        else
            M = M*1000
        end
            
        println("\n O valor do M é: ", M)

        #__________________________________________________________________________________
        # Criar vetor com a coluna b

        b_valores = []

        for i = 1:n_linhas
            push!(b_valores, (tableau_inicial[i,n_colunas]))
        end

        #__________________________________________________________________________________
        # O tableua inicial guarda a matriz até a penultima coluna

        #"minimização" ou "maximização"
        
        tableau_parcial = tableau_inicial[1:n_linhas, 1:n_colunas-1]

        if tipo_otimização == "minimização"
            tableau_parcial[1,:]= tableau_inicial[1, 1:n_colunas-1] * (-1)
        else
            tableau_parcial[1,:] = tableau_inicial[1, 1:n_colunas-1]
        end

        #__________________________________________________________________________________
        # Criar vetor que será contatenado a linha da função objetivo. Esse vetor no fim é transposto

        vetor_Lzero = Float64[]
            
        for i = 2:(n_linhas)
            if vetor_sinal[i] == ">="
                push!( vetor_Lzero, 0.0)
                push!( vetor_Lzero, M)
            elseif vetor_sinal[i] == "="
                push!( vetor_Lzero, M)
            else
                push!( vetor_Lzero, 0.0)
            end
        end

            vetor_Lzero_Transposto = vetor_Lzero'

        ### println(vetor_Lzero_Transposto)
        
        #__________________________________________________________________________________ 
        # o t vai servir para contar a quantidade de colunas da matriz com as folgas ou variáveis artificiais
        #Começa em 2 porque ignora o sinal da FO
        
        t = 0
        
        for i = 2:(n_linhas)
            global t
            if vetor_sinal[i] == ">="
                t = t + 2
            else
                t = t + 1
            end
        end

        #__________________________________________________________________________________ 
        # A nova matriz será criada para adicionar as variáveis artificiais e de folga nas restrições
        #Começa em 2 porque ignora o sinal da FO
        
        novamatriz = zeros(n_linhas-1, t)

        j = 0
        for i = 2:(n_linhas)
            global j
            if vetor_sinal[i] == ">="
                j = j + 1
                novamatriz[i-1,j] = -1
                j = j + 1
                novamatriz[i-1,j] = 1
            else
                j = j + 1
                novamatriz[i-1,j] = 1
            end
        end
        
        #__________________________________________________________________________________ 
        # matriz_Soma será criada para permitir a eliminação das variáveis artificiais da FO

        matriz_Lzero_Novamatriz = [vetor_Lzero_Transposto ;  novamatriz]
        matriz_SuporteSoma = [tableau_parcial matriz_Lzero_Novamatriz b_valores]

        tableau = [tableau_parcial matriz_Lzero_Novamatriz b_valores]

        #println("Isso é a matriz união linha0 e nova matriz: ", matriz_Lzero_Novamatriz)
        #println("Essa matriz que gerará o vetor da FO: ", matriz_SuporteSoma)

        #__________________________________________________________________________________ 
        # Será gerado a nova Linha da FO
        
        (linhas_suporte_Soma, colunas_suporte_Soma) =size(matriz_SuporteSoma)
            
        for i = 2 : linhas_suporte_Soma
            if vetor_sinal[i] == "<="
                matriz_SuporteSoma[i,:] = matriz_SuporteSoma[i,:] * 0
            else
                matriz_SuporteSoma[i,:] = matriz_SuporteSoma[i,:] * (-M)
                end
        end

        vetor_soma_Linhazero = sum(matriz_SuporteSoma, dims=1)

        #println("A linha zero inicial é: ", vetor_soma_Linhazero)
    
        #__________________________________________________________________________________ 
        # Gerar tableau!

        tableau[1,:] = vetor_soma_Linhazero
            
        #println("\n O tableau do problema é: \n")
            
        (line_tableau, colun_tableau)  = size(tableau)

        header_final = String[]

        for j = 1:(colun_tableau)
            push!(header_final, "x$j")
        end

        header_final[colun_tableau, 1] = "b"

        println("\nA matriz com as variáveis de folga/artificial: ")
        pretty_table(tableau, header = header_final, border_crayon = crayon"yellow",formatters = ft_round(2))
  
        ###############################################################################################
        #INICIO DO SIMPLEX
        ###############################################################################################

        (line_tableau, colun_tableau)  = size(tableau)

        iteração = 0

    while minimum(tableau[1,1:(colun_tableau-1)])<0;

        global iteração
        iteração = iteração + 1

        println("\nA iteração é: ",iteração)

         #__________________________________________________________________________________
        #Identificando o Indice da Coluna Pivo 
                
        for j = 1:(colun_tableau-1)
            global indice_coluna
            if minimum(tableau[1,1:(colun_tableau-1)]) == tableau[1, j]
                indice_coluna = j
            end
        end
    
        println("\nO indice da coluna pivô é: ", indice_coluna)
        
        #__________________________________________________________________________________
        #Criando array para o teste da razão minima
        
        a = []

        for  i = 1:line_tableau
            if (tableau[i , colun_tableau] /tableau[i , indice_coluna]<=0)
                push!(a, 1000)
            else    
                push!(a, (tableau[i , colun_tableau] /tableau[i , indice_coluna]))
            end
        end

        ### println("O array com os valores da divisão é:", a)
        
        #__________________________________________________________________________________
        #Achando a linha pivo   
        
        for i = 1:line_tableau
            global indice_linha
            if minimum(a[2:line_tableau,1]) == a[i, 1]
                indice_linha = i
                break
            end
        end
        
        println("O indice da linha pivô é: ", indice_linha)  
        #__________________________________________________________________________________
        #Achando o elemento pivo
        
        elemento_pivo = tableau[indice_linha, indice_coluna]
        
        println("Meu elemento pivo é: ", round(elemento_pivo, digits=2))
        println("\n")
        
        #__________________________________________________________________________________
        #Definição da nova linha pivo
        
        novalinha_pivo = tableau[indice_linha, :] / elemento_pivo

        ### println("A nova linha pivo: ",)
        ### writedlm(stdout, novalinha_pivo)

        #__________________________________________________________________________________
        #Definição das demais linhas
        
        for i = 1:line_tableau
            if i == indice_linha
                tableau[i,:] = novalinha_pivo
            else     
                tableau[i,:] = novalinha_pivo*(-1)*tableau[i,indice_coluna]+tableau[i,:]
            end
        end

            pretty_table(tableau, header = header_final, border_crayon = crayon"yellow",formatters = ft_round(2))
    end
        #__________________________________________________________________________________
        #Sumário Final
        
        println("\nO número de iterações para identificar o valor ótimo foi: ", iteração)
           
        if tipo_otimização == "minimização"
            println("O valor ótimo da função Z é : ",(-1)* round(tableau[1, colun_tableau], digits=2))                
            println("\n")
        else 
            println("O valor ótimo da função Z é : ",round(tableau[1, colun_tableau], digits=2))
            println("\n")
        end

            
        ###############################################################################################
        #DUAL
        ###############################################################################################
        #__________________________________________________________________________________
        #tableau_dual_inicio é o recorte que desconsidera a FO e o b
        
        tableau_dual_inicio = (tableau_inicial[2:n_linhas,1:n_colunas-1])'

        ### println(tableau_dual_inicio)
        #__________________________________________________________________________________
        #Fo_Dual é a FO do dual        
        
        Fo_Dual = zeros(n_linhas-1, 1)

        Fo_Dual = (tableau_inicial[2:n_linhas, n_colunas])'

        #__________________________________________________________________________________
        #b_Dual é o "c" correspondente ao "b" no primal        

        b_Dual_parcial = zeros(1, n_colunas)

        for j = 2:n_colunas
            b_Dual_parcial[1,j] = (tableau_inicial[1,j-1]*(-1))
        end
            
        b_Dual = b_Dual_parcial'

        ### println(b_Dual)

        tableau_dual_parcial = [Fo_Dual;tableau_dual_inicio]

        #__________________________________________________________________________________
        #tableau_dual é o tableau dual do primal      

        tableau_dual = [tableau_dual_parcial b_Dual]
            
        tableau_dual[1,:] = tableau_dual[1,:] *-1

        #__________________________________________________________________________________
        #Ajustando impressão do tableau_dual      
        
        (linhas_dual,colunas_dual) = size(tableau_dual)
    
        header_dual = []

        for j = 1:(colunas_dual)
            push!(header_dual, "y$j")
        end

        header_dual[colunas_dual, 1] = "c"

        println("\nA matriz dual: ")

        pretty_table(tableau_dual, header = header_dual, border_crayon = crayon"yellow",formatters = ft_round(2))

        #__________________________________________________________________________________
        #Ajuste do vetor sinal, para mudança dos sinais de desigualdade  
        
        
        tipo_otimização_dual = []
        vetor_sinal_dual = []
            
        if tipo_otimização == "maximização"
            tipo_otimização_dual = "minimização"
        else
            tipo_otimização_dual = "maximização"
        end
    
        println("\nO tipo de otimização do tableau do dual é de: ",tipo_otimização_dual)
        println("\n")

        for i = 2:n_linhas
            if tipo_otimização == "maximização"
                if vetor_sinal[i,1] == ">="
                    push!(vetor_sinal_dual, "<=")
                elseif vetor_sinal[i,1] == "<="
                    push!(vetor_sinal_dual, ">=")
                else
                    push!(vetor_sinal_dual, "livre")
                end
            else 
                if vetor_sinal[i,1] == ">="
                    push!(vetor_sinal_dual, ">=")
                elseif vetor_sinal[i,1] == "<="
                    push!(vetor_sinal_dual, "<=")
                else
                    push!(vetor_sinal_dual, "livre")
                end 
            end 
        end

        #pretty_table(vetor_sinal_dual)

            
        Solucao_Dual_inicial = []
        Solucao_Dual_inicial = tableau[1, (n_colunas:colun_tableau-1)]

        #println("Solução inicial", Solucao_Dual_inicial)

        Solucao_Dual_parcial = hcat(Solucao_Dual_inicial...)

        #println("Nova Solução inicial", Solucao_Dual_parcial)

        #println(vetor_Lzero_Transposto)

            #__________________________________________________________________________________
        # Identificar as colunas que são variável de folga e desconsiderar as artificiais


        Solucao_Dual_final = []
            
        for j = (1:length(vetor_Lzero_Transposto))
            if vetor_Lzero_Transposto[1,j] == 0
                push!(Solucao_Dual_final, Solucao_Dual_parcial[1, j])  
            end
        end

        println("A solução Dual é: \n")

        variavel_dual = []

        for i = 1:(n_linhas-1)
            push!(variavel_dual, "y$i")
        end

        #println(variavel_dual)

        pretty_table(hcat(Solucao_Dual_final...)', row_names = variavel_dual, header = ["Valor"], border_crayon = crayon"yellow",formatters = ft_round(2))

         

        ###############################################################################################
        #ANALISE DE SENSIBLIDADE
        ###############################################################################################

        #Baseado na equação da b = S*DeltaB+b* >=0    
    
        #Cortar o tableau Final

        S_inicial = []

        S_inicial = tableau[2:n_linhas, (n_colunas:colun_tableau-1)]

        ## vetor_Lzero_Transposto aponta caso diferente de zero, a posição da coluna da variável artificial

        S_parcial = []
            
        for j = (1:length(vetor_Lzero_Transposto))
            if vetor_Lzero_Transposto[1,j] == 0
                push!(S_parcial, S_inicial[:, j])  
            end
        end

        #println(S_parcial)

        S_fim =[]
        S_fim  = hcat(S_parcial...)

        (lin_S_fim, col_S_fim) =size(S_fim)

        DeltaB = zeros(col_S_fim,col_S_fim) + I
            
        #pretty_table(DeltaB)
        #pretty_table(S_fim)

        #O b_sensibilidade é igual ao b da matriz sem variáveis de folga e artificial
        #colun_tableau é o número de colunas do tableau
            
        b_sensibilidade =[]
        b_sensibilidade = tableau[2:n_linhas,colun_tableau]

        #sinais_inicio da suporte para equação realizada!

        sinais_inicio = Array{String}(undef, col_S_fim, col_S_fim)

        for j =(1:col_S_fim)
            for i = (1:col_S_fim)
                sinais_inicio[i,j] =">="
            end
        end

        #println("Sinais inicio", sinais_inicio)

        sinais_fim = sinais_inicio

        for i = 1:col_S_fim
            for j =1:col_S_fim
                if S_fim[i,j] < 0
                    sinais_fim[i, j] = "<="
                else
                    sinais_fim[i,j] = ">="
                end
            end
        end
            
        #println("sinais fim", sinais_fim)

        b_sensibilidade = b_sensibilidade*-1
        
        Matriz_parcial = []

        #println("S_Fim", S_fim)
        #println("DeltaB", DeltaB)
        
        Matriz_parcial= S_fim * DeltaB

        #println("matriz parcial", Matriz_parcial)
        
        Matriz_fim =zeros(col_S_fim,col_S_fim)
        
        for i = 1:col_S_fim
            for j =1:col_S_fim
                if Matriz_parcial[i,j] < 0 
                    Matriz_fim[i,j] = b_sensibilidade[i,1]/Matriz_parcial[i,j]
                else
                    Matriz_fim[i,j] = b_sensibilidade[i,1]/Matriz_parcial[i,j]
                end
            end
        end
            
            
        #println("matriz fim", Matriz_fim)
            
        vetor_maior_igual = zeros(col_S_fim,col_S_fim)
        vetor_menor_igual = zeros(col_S_fim,col_S_fim)
 
        for j = 1:col_S_fim
            for i = 1:col_S_fim
                if sinais_fim[i,j] == ">="
                    vetor_maior_igual[i,j] = Matriz_fim[i,j]
                else
                    vetor_maior_igual[i,j] = -Inf
                end     
            end 
        end

        #println("vetor maior igual", vetor_maior_igual)

        for j = 1:col_S_fim
            for i = 1:col_S_fim
                if sinais_fim[i,j] == "<="
                    vetor_menor_igual[i,j] = Matriz_fim[i,j]
                else
                    vetor_menor_igual[i,j]  = Inf
                end
            end 
        end
            
            #println("vetor menor igual", vetor_menor_igual)
            
        range_inferior = []
        range_superior = []
        
        for j = 1:col_S_fim
            push!(range_inferior, maximum(vetor_maior_igual[:,j]))
            push!(range_superior, minimum(vetor_menor_igual[:,j]))
        end


        #println("O variação do DeltaB no limite inferior pode ser: ", range_inferior)
            
        #println("O variação do DeltaB no limite superior pode ser: ", range_superior)


        teste = [range_superior, range_inferior]

        teste_fim  = hcat(teste...)

        data =[]

        for i =1:col_S_fim
            push!(data, "Restrição $i")
        end

        data_fim = (teste_fim)

        println("\nAnálise de sensibilidade nas restrições: \n")

        println("\nAs variações do lado direito das restrições podem ser: \n")

        pretty_table(data_fim, row_names = data, header = ["Incremento permitido", "Decremento permitido"], border_crayon = crayon"yellow",formatters = ft_round(2))


