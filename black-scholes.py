# -*- coding: utf-8 -*-
""" 
Augusto Vaccarelli Costa - NUSP 10770197
Italo Roberto Lui - NUSP 12553991
MAP3122 : A2 - Modelo de Black-Scholes

# O programa abaixo aplica o Método das Linhas, Euler Implícito e Interpolação Polinomial (Quadrática),
# à um problema representado pelo modelo de precificação de ativos Black-Scholes, dados os parâmetros σ, K, T e r.
# O resultado é a interpolação quadrática do preço da opção, exibida em um gráfico.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from statistics import NormalDist

################################### Variáveis e valores exatos no modelo ##########################################
def tau(T, t):
    # Mudança de variável
    return (T - t)

def S(sigma, K, r, x, tau):
    # Preço da ação/ativo em determinado tempo
    return (K * math.e**(x - ((r - ((sigma**2)/2)) * tau)))

def dezinho_1(sigma, x, tau):
    return ((x + (sigma**2 * tau))/(sigma * math.sqrt(tau)))

def dezinho_2(sigma, x, tau):
    return (x/(sigma * math.sqrt(tau)))

def dezao_1(sigma, K, r, x, tau):
    return ((math.log(S(sigma, K, r, x, tau)/K) + ((r + (sigma**2)/2) * tau))/(sigma * math.sqrt(tau)))

def dezao_2(sigma, K, r, x, tau):
    return (dezao_1(sigma, K, r, x, tau) - (sigma * math.sqrt(tau)))

def u(sigma, K, x, tau):
    # Mudança de variável
    return (K * math.e**(x + (sigma**2 * tau/2)) * NormalDist().cdf(dezinho_1(sigma, x, tau)) - (K * NormalDist().cdf(dezinho_2(sigma, x, tau))))

def V(sigma, K, r, x, tau):
    # Preço da opção de compra europeia
    return (S(sigma, K, r, x, tau) * NormalDist().cdf(dezao_1(sigma, K, r, x, tau)) - (K * math.e**((-1 * r) * tau) * NormalDist().cdf(dezao_2(sigma, K, r, x, tau))))

def difdiv(j, *tabela):
    # Diferença dividida de ordem j da função tabelada, sendo que tabela[0] tem os valores de x e tabela[1] os valores de y.
    numerador = 0
    denominador = 0
    i = 0

    for i in range(0, j + 1):
        pass 

def metodoLinhas(sigma, K, r, T, L, N, M):
    # Aplica o Método das Linhas à equação de Black-Scholes no formato da equação de difusão do calor.
    # Cada linha da matriz representa uma posição fixa em x, e cada coluna representa um momento tau.
    deltaX   = L/N
    deltaTau = T/M
    i = 0
    j = 0
    u = np.zeros((N + 1, M + 1))

    for i in range (0, N + 1):
        for j in range (0, M + 1):
            if (j == 0):
                if ((K * math.e**((-1 * L) + (i * deltaX)) - 1) < 0):
                    u[i][j] = 0
                else:
                    u[i][j] = (K * math.e**((-1 * L) + (i * deltaX)) - 1)
            elif (i == 0):
                u[i][j] = 0
            elif (i == N):
                u[i][j] = (K * math.e**(L + ((sigma**2 * tau(T, j * deltaTau))/2)))
            else:
                u[i][j] = (u[i][j-1] + ((deltaTau/(deltaX**2)) * (sigma**2/2) * (u[i-1][j-1] - (2 * u[i][j-1]) + u[i+1][j-1])))

    return u

def main():
    # Dados Iniciais do Modelo
    sigma = 0 # Volatilidade (adim.)
    K     = 0 # Preço de exercício da opção (u.m.)
    T     = 0 # Tempo de vencimento / data de exercício da opção (anos)
    r     = 0 # Taxa de juros livre de risco anualizada

    # Parâmetros do programa
    L = 10    # Subdomínio de x usado para o Método das Linhas é [-L, L]
    N = 1000  # Número de passos de integração de x para o Método das Linhas
    M = 100   # Número de passos de integração de tau para o Método das Linhas

    print("Insira os dados necessários para a precificação da opção.\n")
    print("\n")
    sigma = float(input("Volatilidade do ativo: "))
    K     = float(input("Preço de exercício da opção: "))
    T     = float(input("Data de exerício da opção (em anos): "))
    r     = float(input("Taxa de juros livre de risco (anualizada): "))

    for M in (100, 200, 300, 400, 500):
        u = metodoLinhas(sigma, K, r, T, L, N, M)

main()