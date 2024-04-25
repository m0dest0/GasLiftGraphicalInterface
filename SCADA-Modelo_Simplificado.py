# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 11:44:57 2024

@author: Pedro
"""
from pyModbusTCP.client import ModbusClient
from pyModbusTCP import utils
import numpy as np
from scipy.integrate import odeint 
import time

class FloatModbusClient(ModbusClient):
    def read_float(self, address, number=1):
        reg_l = self.read_holding_registers(address, number * 2)
        if reg_l:
            return [utils.decode_ieee(f) for f in utils.word_list_to_long(reg_l)]
        else:
            return None

    def write_float(self, address, floats_list):
        b32_l = [utils.encode_ieee(f) for f in floats_list]
        b16_l = utils.long_list_to_word(b32_l)
        return self.write_multiple_registers(address, b16_l)

def model(x,t,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17):  
    
    x = x.tolist()
    
    # Constantes
    
    Ta = T1  # Temperatura da região anular [K]
    Tt = T2  # Temperatura no tubo [K]
    ro = T3  # Densidade do óleo no reservatório [kg/m3]
    Pr = T4 * 1000  # Pressao no reservatório longe da cabeça do poço [Pa];
    Ps = T5 * 1000  # Pressão no manifold [Pa];
    wgc = T6 # Vazão do gás que entra na região anular [m³/s]
    Cpc = 1.655e-3 * T7/100 # [m2]
    Civ = 1.5e-4 * T8/100# [m2]
    Cr = 2.623e-4 * T9/100 # [m2]
    La = T10  # Comprimento da região anular [m]
    Va = La * np.pi * ((T15 - T11)**2)/4   # Volume da região anular [m3]
    Lt = T12  # Comprimento do tubo [m]
    Lr = T13  # Distância do reservatório até o ponto de injeção [m]
    Ar = np.pi * (T14**2)/4   # Área da seção transversal abaixo do ponto injeção [m2]
    At = np.pi * (T15**2)/4    # Área da seção transversal acima do ponto injeção [m2]
    M = T16  # Massa molar do gás [kg/mol]
    Vo = 1/ro  # volume específico do óleo
    
    g = 9.81  # m/s2
    R = 8.314  # J/mol.K
    
    # Pressões em diferentes pontos

    Pai=((R*Ta/(Va*M) + g * La/Va )* x[0])  # Pressão no anular no ponto de injeção
    Pt = (R * Tt /M ) * (x[1] / (Lt * At - Vo * x[2]))  # Pressão no topo do tubo
    Pti = Pt + ( g/At )* (x[1] + x[2])   # Pressão no tubo no ponto de injeção
    Ptb = Pti + ro * g * Lr  # Pressão no fundo do tubo
    
    # Densidades
    ro_ai = (M * Pai) / (R * Ta)  # Densidade do gás na regiao anular no ponto de injeção
    ro_m = ( x[1] + x[2] ) / (Lt * Ar)  # densidade da mistura óleo/gás no tubo

    # Vazões

    m1 = max(0, Pai - Pti)
    m2 = max(0, Pr - Ptb)
    m3 = max(0,Pt - Ps )

    wiv = Civ * np.sqrt(ro_ai * m1)  # Vazão de gás na Válvula de gas lift
    wpc = Cpc * np.sqrt(ro_m * m3)  # Vazão na válvula Choke (Vazão total na cabeça do poço)
    wpg = x[1] * wpc / (x[1] + x[2])  # Vazão de gás na cabeça do poço
    wpo = x[2] * wpc / (x[1] + x[2])  # Vazão de óleo na cabeça do poço
    wr = Cr * np.sqrt(ro * m2)  # Vazão de óleo na fronteira reservatório / tubo

    # Sistema de EDO'S

    dxdt = np.zeros(3) # Cria uma matriz com cada equação de balanço na forma diferencial

    dxdt[0] = wgc - wiv
    dxdt[1] = wiv - wpg
    dxdt[2] = wr - wpo

    return dxdt

def plot(x,t,I): 
    
    x = x.tolist()
    
    # Constantes
    
    Ta = I[1]  # Temperatura da região anular [K]
    Tt = I[2]  # Temperatura no tubo [K]
    ro = I[3] # Densidade do óleo no reservatório [kg/m3]
    Pr = I[4] * 1000  # Pressao no reservatório longe da cabeça do poço [Pa];
    Ps = I[5] * 1000 # Pressão no manifold [Pa];
    Cpc = 1.655e-3 * I[7]/100 # [m2]
    Civ = 1.5e-4 * I[8]/100 # [m2]
    Cr = 2.623e-4 * I[9]/100
    La = I[10]  # Comprimento da região anular [m]
    Va = La * np.pi * ((I[15] - I[11])**2)/4   # Volume da região anular [m3]
    Lt = I[12]  # Comprimento do tubo [m]
    Lr = I[13]  # Distância do reservatório até o ponto de injeção [m]
    Ar = np.pi * (I[14]**2)/4   # Área da seção transversal abaixo do ponto injeção [m2]
    At = np.pi * (I[15]**2)/4    # Área da seção transversal acima do ponto injeção [m2]
    M = I[16]  # Massa molar do gás [kg/mol]
    Vo = 1/ro  # volume específico do óleo
    
    g = 9.81  # m/s2
    R = 8.314  # J/mol.K
   
    # Pressões em diferentes pontos

    Pai=((R*Ta/(Va*M) + g * La/Va )* x[0])  # Pressão no anular no ponto de injeção
    Pt = (R * Tt /M ) * (x[1] / (Lt * At - Vo * x[2]))  # Pressão no topo do tubo
    Pti = Pt + ( g/ At )* (x[1] + x[2])   # Pressão no tubo no ponto de injeção
    Ptb = Pti + ro * g * Lr  # Pressão no fundo do tubo
    
    # Densidades
    ro_ai = (M * Pai) / (R * Ta)  # Densidade do gás na regiao anular no ponto de injeção
    ro_m = ( x[1] + x[2] ) / (Lt * Ar)  # densidade da mistura óleo/gás no tubo

    # Vazões

    m1 = max(0, Pai - Pti)
    m2 = max(0, Pr - Ptb)
    m3 = max(0,Pt - Ps )

    wiv = Civ * np.sqrt(ro_ai * m1)  # Vazão de gás na Válvula de gas lift
    wpc = Cpc * np.sqrt(ro_m * m3)  # Vazão na válvula Choke (Vazão total na cabeça do poço)
    wpg = x[1] * wpc / (x[1] + x[2])  # Vazão de gás na cabeça do poço
    wpo = x[2] * wpc / (x[1] + x[2])  # Vazão de óleo na cabeça do poço
    wr = Cr * np.sqrt(ro * m2)  # Vazão de óleo na fronteira reservatório / tubo
    
    return np.array([wr, wpc, wpg, wpo, Pt, Pti, Ptb, Pai, ro_ai, ro_m, wiv])

c = FloatModbusClient(host='localhost', port=5022, auto_open=True) # Localização do Servidor 
  
while True: # Loop de funcionamento do programa

    # Condições iniciais
    x0 = [4350.1, 10951, 86038]
    t0 = 0
    delta_t = 1 
    
    B = False # Variável Booleana que serve como chave da simulação
    
    b = c.read_float(0,1) # Leitura do valor do botão ON/OFF no ScadaBR
    
    if b[0] < 0: # Se o botão estiver ligado, a variável muda de estado e inicia-se a simulação
        B = True
        
    while B == True: # Loop de simulação
          
        I = c.read_float(0,18) # Leitura dos inputs no ScadaBR
        
        I = tuple(I) # Conversão dos valores lidos para manipulação matemática
        
        t = np.linspace(t0, t0 + delta_t) # Array de tempo
        
        sol = odeint( model, x0, t, I) # Integração do Modelo em t
        X,Y,Z = sol. transpose() # Transposição de uma matriz 3 x N para N x 3
        x0 = [X[-1],Y[-1],Z[-1]] # Atualização da condição inicial com os últimmos valores
        
        V = plot(sol[-1],t,I) # Função para plotar as variáveis dentro do modelo

        
        t0 = t0 + delta_t # Atualização da condição inicial
        
        c.write_float(0,[I[0],I[1],I[2],I[3],I[4],I[5],I[6],I[7],I[8],I[9],I[10],
                         I[11],I[12],I[13],I[14],I[15],I[16],I[17],V[0],
                         V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],t0]) #Escrita no ScadaBR dos inputs seguido dos outputs
                
        # Verificando se o botão ON/OFF foi desligado
        
        if I[0] > 0:  # Se sim,verifica as condições da variável manter histórico
            
            if I[17] > 0: # Se ela estiver desligada, sai-se do loop de simulação e zeram-se as condições iniciais e de tempo
                B = False
            else: # Se estiver ligada, entra-se nesse loop para manter as condições iniciais até que seja ligada novamente a simulação
                while I[17] < 0 and I[0] > 0:
                    I = c.read_float(0,19)
                    
                    time.sleep(1) # Limitação do tempo de execução
                     
        time.sleep(1) # Limitação do tempo de execução
    
    time.sleep(1) # Limitação do tempo de execução
    
           
        
      
    


