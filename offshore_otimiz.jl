# Offshore_otimiz

#Trabalho Final - Otimização Linear

using JuMP 
using Cbc
using PrettyTables

    n_vertices = 7       #É a soma do nó do Porto e das N Turbinas
    #n_turbinas

    #p é prêmio da visita ao local
    
    penalidade = [0 200 200 200 200 1000 200]

    premio = [0 200 200 200 200 200 200]

    Tempo_Manutenção = [0 1 1 1 1 3 1]
    
    Tempo_Total = 8

    pmin = 0

    ## Definição do premio minimo

    for i = 1: length(premio)
        global pmin
        pmin = pmin + premio[i]
    end

    pmin = pmin * 1/3

    println(pmin)

    ## Custo de movimentação

    custo = Float64[]

    custo = [0 50 100 150 200 250 300; 50 0 50 100 150 200 250; 100 50 0 50 100 150 200; 150 100 50 0 50 100 150; 200 150 100 50 0 50 100; 250 200 150 100 50 0 50; 300 250 200 150 100 50 0]
 
    ## Tempo de deslocamento entre nós
    
    Tempo_Deslocamento = Float64[]
    
    Tempo_Deslocamento = [0 0.25 0.5 0.75 1 1.25 1.5; 0.25 0 0.25 0.5 0.75 1 1.25; 0.5 0.25 0 0.25 0.5 0.75 1; 0.75 0.5 0.25 0 0.25 0.5 0.75; 1 0.75 0.5 0.25 0 0.25 0.5; 1.25 1 0.75 0.5 0.25 0 0.25; 1.5 1.25 1 0.75 0.5 0.25 0]

    model = Model(Cbc.Optimizer)

    @variable(model, x[i=1:n_vertices,j=1:n_vertices], Bin)
    @variable(model, y[i=1:n_vertices], Bin)

    @variable(model, f[i=1:n_vertices, j=1:n_vertices]>=0)
    #Nao deveria existir. Ajustar! Não precisa ser inteira

    ## Função objetivo 

    Custo_Total      = sum(custo[i,j] * x[i,j] for i=1:n_vertices, j=1:n_vertices if j!=n_vertices)
    Penalidade_Total = sum(penalidade[i]*(1-y[i]) for i=1:n_vertices)

        @objective(model, Min, Custo_Total + Penalidade_Total)

    ## Restrição (1)
    for i = 1:n_vertices
        @constraint(model, sum(x[i,j] for j=1:n_vertices if j != i)==y[i])
    end

    ## Restrição (2)
    for j = 1:n_vertices
        @constraint(model, sum(x[i,j] for i=1:n_vertices if i != j)==y[j])
    end

    ## Restrição (3)
        Premio = @constraint(model, sum(premio[i]*y[i] for i=1:n_vertices) >= pmin)

    ## Restrição (4)
    for i = 2:n_vertices
        @constraint(model, sum(f[j,i] for j=1:n_vertices if j!=i) - sum(f[i,j] for j=1:n_vertices if j!=i ) == y[i])
    end

    ## Restrição (5)
    for i = 1:n_vertices
    for j = 1:n_vertices
    if  j!= i
                @constraint(model, f[i,j] <= (n_vertices-1)*x[i,j])      
            end
        end
    end

    ## Restrição (6) - Adição de Tempo

    Tempo_Deslocamento = sum(Tempo_Deslocamento[i,j]*x[i,j] for i=1:n_vertices, j= 1:n_vertices)
    Tempo_Manutenção   = sum(Tempo_Manutenção[i]*y[i] for i = 1:n_vertices)

    Tempo_Operacao = @constraint(model, Tempo_Deslocamento + Tempo_Manutenção<= Tempo_Total)

    
    ## Restrição adicional = garantir que o nó do Porto será visitado
        @constraint(model, y[1]==1)

    ## Otimização do modelo    

    print(model)

    optimize!(model)

    pretty_table(value.(x))
    pretty_table(value.(y))
    pretty_table(value.(f))

    status       = JuMP.optimize!(model)
    status       = termination_status(model)
    PremioMaximo = JuMP.value(Premio)

    println("STATUS: ", status, " ---------------------------")

    solvalue = objective_value(model)
    println("Função Objetivo ",solvalue)

    println("Premio Máximo: ",PremioMaximo)

    println("Custo: ", JuMP.value(Custo_Total))
    println("Penalidade: ", JuMP.value(Penalidade_Total))
    println("Tempo de Operação: ",JuMP.value(Tempo_Operacao))
    println("\n")

    ## Impressão da rota realizada
    
    println("A rota será detalhada abaixo: \n")

    for i = 1 : n_vertices
        if i == 1
            print("Do Porto      ")
        else
            print("Da Turbina : ", (i-1))
        end  
        for j = 1:n_vertices
            if (value.(x[i,j]) != 0)
                if j == 1
                    print(" para o Porto")
                else
                print(" para a Turbina: ", (j-1))    
                end
            end
        end
        println()
    end
    println("\n")