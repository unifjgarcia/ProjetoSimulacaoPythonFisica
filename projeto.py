import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import messagebox
from matplotlib.animation import FuncAnimation


pi = 3.14
h = 6.626E-34
mE = 9.109E-31
mP = 1.6726e-27
eEV = 1.602e-19
c = 3E8
global frame_grafico
frame_grafico = None

fig, ax = plt.subplots()
ax.set_xlim(0, 10)
ax.set_ylim(-1.5, 1.5)
line, = ax.plot([], [], lw=2)

def densidade_pontoX(x, n, L):
    return np.sqrt(2/L) * np.sin(n * np.pi * x / L)

def calculosProton(resultado_text, limite_inf, limite_sup, n_energia, n2_energia, largura_poco):
    # Conversão de tipos
    limite_inf = float(limite_inf)
    limite_sup = float(limite_sup)
    n_energia = int(n_energia)
    n2_energia = int(n2_energia)
    largura_poco = float(largura_poco)
    
    # Cálculo da função de onda inicial
    p1 = 2 / largura_poco
    ptotal = np.sqrt(p1)
    ntotal = (n_energia * np.pi) / largura_poco
    resultado_formatado = "{:.3e}".format(ptotal)
    resultado_formatado2 = "{:.3e}".format(ntotal)
    resultado_text.insert(tk.END, f"Sua função de onda é {resultado_formatado} sen({resultado_formatado2}) * x\n")

    # Cálculo da função de onda final
    p1 = 2 / largura_poco
    ptotal = np.sqrt(p1)
    n2total = (n2_energia * np.pi) / largura_poco
    resultado_formatado = "{:.3e}".format(ptotal)
    resultado_formatado2 = "{:.3e}".format(n2total)
    resultado_text.insert(tk.END, f"Sua função de onda é {resultado_formatado} sen({resultado_formatado2}) * x\n")

    # Cálculo da energia total inicial do próton
    e = np.multiply(8, mP)
    e2 = np.multiply(e, largura_poco**2)
    e3 = np.divide(h**2, e2)
    etotal1 = np.multiply(e3, n_energia**2)
    etotal2 = np.divide(etotal1, eEV)
    resultado_formatado = "{:.4e}".format(etotal1)
    resultado_formatado2 = "{:.4e}".format(etotal2)
    resultado_text.insert(tk.END, f"Sua energia total é {resultado_formatado}J = {resultado_formatado2}) eV\n")

    # Cálculo da energia total final do próton
    enf = np.multiply(8, mP)
    enf2 = np.multiply(enf, largura_poco**2)
    enf3 = np.divide(h**2, enf2)
    etotalnf1 = np.multiply(enf3, n2_energia**2)
    etotalnf2 = np.divide(etotalnf1, eEV)
    resultado_formatado = "{:.4e}".format(etotalnf1)
    resultado_formatado2 = "{:.4e}".format(etotalnf2)
    resultado_text.insert(tk.END, f"Sua energia total é {resultado_formatado}J = {resultado_formatado2}) eV\n")

    # Cálculo da energia absorvida pelo fóton
    e = 8 * mP
    e2 = e * largura_poco**2
    e3 = h**2 / e2
    etotal1 = e3 * n_energia**2
    etotal2 = etotal1 / eEV
    
    en2 = 8 * mP
    e2n2 = en2 * largura_poco**2
    e3n2 = h**2 / e2n2
    etotaln2 = e3n2 * n2_energia**2
    etotal2n2 = etotaln2 / eEV
    
    energia_absorvida_J = etotaln2 - etotal1
    energia_absorvida_eV = etotal2n2 - etotal2
    
    resultado_formatado = "{:.4e}".format(abs(energia_absorvida_J))
    resultado_formatado2 = "{:.4e}".format(abs(energia_absorvida_eV))
    
    resultado_text.insert(tk.END, f"Sua energia absorvida é {resultado_formatado} J ou {resultado_formatado2} eV\n")

    # Cálculo do comprimento de onda
    conta1 = np.multiply(h,c)
    cdeonda = np.divide(conta1, energia_absorvida_J)
    resultado_formatado = "{:.4e}".format(abs(cdeonda))
    resultado_text.insert(tk.END, f"Seu comprimento de onda é {resultado_formatado} m\n")

    # Cálculo da frequencia
    frequencia = np.divide(c,cdeonda)
    resultado_formatado = "{:.4e}".format(abs(frequencia))
    resultado_text.insert(tk.END, f"Sua frequencia é {resultado_formatado} Hz\n")

    # Cálculo da velocidade inicial
    conta1 = np.multiply(etotal1,2)
    conta2 = np.divide(conta1,mP)
    velocidadeInicial = sqrt(conta2)
    resultado_formatado = "{:.4e}".format(velocidadeInicial)
    resultado_text.insert(tk.END, f"Sua velocidade inicial é {resultado_formatado} m/s\n")

    # Cálculo da velocidade final
    conta1 = np.multiply(etotaln2,2)
    conta2 = np.divide(conta1,mP)
    velocidadeFinal = sqrt(conta2)
    resultado_formatado = "{:.4e}".format(velocidadeFinal)
    resultado_text.insert(tk.END, f"Sua velocidade final é {resultado_formatado} m/s\n")

    # Cálculo De Broglie inicial
    conta1 = np.multiply(velocidadeInicial,mP)
    cdeondaDeBroglie = np.divide(h,conta1)
    resultado_formatado = "{:.4e}".format(cdeondaDeBroglie)
    resultado_text.insert(tk.END, f"Seu comprimento de onda inicial de De Broglie é {resultado_formatado} m/s\n")

    # Cálculo De Broglie final
    conta1 = np.multiply(velocidadeFinal,mP)
    cdeondaDeBroglie = np.divide(h,conta1)
    resultado_formatado = "{:.4e}".format(cdeondaDeBroglie)
    resultado_text.insert(tk.END, f"Seu comprimento de onda inicial de De Broglie é {resultado_formatado} m/s\n")
    # Cálculo da probabilidade
    integral = 0
    precisao_integral = 1000
    dx = (limite_sup - limite_inf) / precisao_integral  
    
    for i in range(precisao_integral):
        x = limite_inf + i * dx
        integral += np.abs(densidade_pontoX(x, n_energia, largura_poco))**2 * dx

    porcentagem = integral * 100
    resultado_text.insert(tk.END, f"A probabilidade é {'{:.2e}'.format(abs(porcentagem))} %\n")

    integral = 0
    precisao_integral = 1000
    dx = (limite_sup - limite_inf) / precisao_integral

    for i in range(precisao_integral):
        x = limite_inf + i * dx
        integral += np.abs(densidade_pontoX(x, n2_energia, largura_poco))**2 * dx

    porcentagem2 = integral * 100
    resultado_text.insert(tk.END, f"A probabilidade é {'{:.2e}'.format(abs(porcentagem2))} %\n")

def calculosEletron(resultado_text, limite_inf, limite_sup, n_energia, n2_energia, largura_poco):
    # Conversão de tipos
    limite_inf = float(limite_inf)
    limite_sup = float(limite_sup)
    n_energia = int(n_energia)
    n2_energia = int(n2_energia)
    largura_poco = float(largura_poco)
    
    # Cálculo da função de onda inicial
    p1 = 2 / largura_poco
    ptotal = np.sqrt(p1)
    ntotal = (n_energia * np.pi) / largura_poco
    resultado_formatado = "{:.3e}".format(ptotal)
    resultado_formatado2 = "{:.3e}".format(ntotal)
    resultado_text.insert(tk.END, f"Sua função de onda é {resultado_formatado} sen({resultado_formatado2}) * x\n")

    # Cálculo da função de onda final
    p1 = 2 / largura_poco
    ptotal = np.sqrt(p1)
    n2total = (n2_energia * np.pi) / largura_poco
    resultado_formatado = "{:.3e}".format(ptotal)
    resultado_formatado2 = "{:.3e}".format(n2total)
    resultado_text.insert(tk.END, f"Sua função de onda é {resultado_formatado} sen({resultado_formatado2}) * x\n")

    # Cálculo da energia total inicial do eletron
    e = np.multiply(8, mE)
    e2 = np.multiply(e, largura_poco**2)
    e3 = np.divide(h**2, e2)
    etotal1 = np.multiply(e3, n_energia**2)
    etotal2 = np.divide(etotal1, eEV)
    resultado_formatado = "{:.4e}".format(etotal1)
    resultado_formatado2 = "{:.4e}".format(etotal2)
    resultado_text.insert(tk.END, f"Sua energia total é {resultado_formatado}J = {resultado_formatado2}) eV\n")

    # Cálculo da energia total final do eletron
    enf = np.multiply(8, mE)
    enf2 = np.multiply(enf, largura_poco**2)
    enf3 = np.divide(h**2, enf2)
    etotalnf1 = np.multiply(enf3, n2_energia**2)
    etotalnf2 = np.divide(etotalnf1, eEV)
    resultado_formatado = "{:.4e}".format(etotalnf1)
    resultado_formatado2 = "{:.4e}".format(etotalnf2)
    resultado_text.insert(tk.END, f"Sua energia total é {resultado_formatado}J = {resultado_formatado2}) eV\n")

    # Cálculo da energia absorvida pelo fóton
    e = 8 * mE
    e2 = e * largura_poco**2
    e3 = h**2 / e2
    etotal1 = e3 * n_energia**2
    etotal2 = etotal1 / eEV
    
    en2 = 8 * mE
    e2n2 = en2 * largura_poco**2
    e3n2 = h**2 / e2n2
    etotaln2 = e3n2 * n2_energia**2
    etotal2n2 = etotaln2 / eEV
    
    energia_absorvida_J = etotaln2 - etotal1
    energia_absorvida_eV = etotal2n2 - etotal2
    
    resultado_formatado = "{:.4e}".format(abs(energia_absorvida_J))
    resultado_formatado2 = "{:.4e}".format(abs(energia_absorvida_eV))
    
    resultado_text.insert(tk.END, f"Sua energia absorvida é {resultado_formatado} J ou {resultado_formatado2} eV\n")

    # Cálculo do comprimento de onda
    conta1 = np.multiply(h,c)
    cdeonda = np.divide(conta1, energia_absorvida_J)
    resultado_formatado = "{:.4e}".format(abs(cdeonda))
    resultado_text.insert(tk.END, f"Seu comprimento de onda é {resultado_formatado} m\n")

    # Cálculo da frequencia
    frequencia = np.divide(c,cdeonda)
    resultado_formatado = "{:.4e}".format(abs(frequencia))
    resultado_text.insert(tk.END, f"Sua frequencia é {resultado_formatado} Hz\n")

    # Cálculo da velocidade inicial
    conta1 = np.multiply(etotal1,2)
    conta2 = np.divide(conta1,mE)
    velocidadeInicial = sqrt(conta2)
    resultado_formatado = "{:.4e}".format(velocidadeInicial)
    resultado_text.insert(tk.END, f"Sua velocidade inicial é {resultado_formatado} m/s\n")

    # Cálculo da velocidade final
    conta1 = np.multiply(etotaln2,2)
    conta2 = np.divide(conta1,mE)
    velocidadeFinal = sqrt(conta2)
    resultado_formatado = "{:.4e}".format(velocidadeFinal)
    resultado_text.insert(tk.END, f"Sua velocidade final é {resultado_formatado} m/s\n")

    # Cálculo De Broglie inicial
    conta1 = np.multiply(velocidadeInicial,mE)
    cdeondaDeBroglie = np.divide(h,conta1)
    resultado_formatado = "{:.4e}".format(cdeondaDeBroglie)
    resultado_text.insert(tk.END, f"Seu comprimento de onda inicial de De Broglie é {resultado_formatado} m/s\n")

    # Cálculo da probabilidade
    integral = 0
    precisao_integral = 1000
    dx = (limite_sup - limite_inf) / precisao_integral  
    
    for i in range(precisao_integral):
        x = limite_inf + i * dx
        integral += np.abs(densidade_pontoX(x, n_energia, largura_poco))**2 * dx

    porcentagem = integral * 100
    resultado_text.insert(tk.END, f"A probabilidade é {'{:.2e}'.format(abs(porcentagem))} %\n")


def plot_grafico(limite_inf, limite_sup, n_energia, n2_energia, largura_poco):
    # Conversão de tipos
    limite_inf = float(limite_inf)
    limite_sup = float(limite_sup)
    n_energia = int(n_energia)
    n2_energia = int(n2_energia)
    largura_poco = float(largura_poco)
    
    # Gerar dados para os gráficos
    x = np.linspace(limite_inf, limite_sup, 1000)
    y_inicial = np.sin(n_energia * np.pi / largura_poco * x)
    y_final = np.sin(n2_energia * np.pi / largura_poco * x)

     # Cálculo da probabilidade para o nível inicial
    integral_inicial = 0
    precisao_integral = 1000
    dx = (limite_sup - limite_inf) / precisao_integral
    
    for i in range(precisao_integral):
        xi = limite_inf + i * dx
        integral_inicial += np.abs(densidade_pontoX(xi, n_energia, largura_poco))**2 * dx

    probabilidade_inicial = integral_inicial * 100

    # Cálculo da probabilidade para o nível final
    integral_final = 0
    
    for i in range(precisao_integral):
        xi = limite_inf + i * dx
        integral_final += np.abs(densidade_pontoX(xi, n2_energia, largura_poco))**2 * dx

    probabilidade_final = integral_final * 100

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
     # Plotar gráfico da função de onda inicial
    axs[0, 0].plot(x, y_inicial, label='Função de Onda Inicial')
    axs[0, 0].set_title('Função de Onda Inicial')
    axs[0, 0].set_xlabel('Posição')
    axs[0, 0].set_ylabel('Função de Onda')
    axs[0, 0].legend()

    # Plotar gráfico da função de onda final
    axs[0, 1].plot(x, y_final, label='Função de Onda Final')
    axs[0, 1].set_title('Função de Onda Final')
    axs[0, 1].set_xlabel('Posição')
    axs[0, 1].set_ylabel('Função de Onda')
    axs[0, 1].legend()

    axs[1, 0].plot(x, y_inicial**2, label='Probabilidade Inicial')
    axs[1, 0].set_title('Probabilidade - Nível Inicial')
    axs[1, 0].set_xlabel('Posição')
    axs[1, 0].set_ylabel('Probabilidade (%)')
    axs[1, 0].legend()

    # Plotar gráfico da probabilidade para o nível final
    axs[1, 1].plot(x, y_final**2, label='Probabilidade Final')
    axs[1, 1].set_title('Probabilidade - Nível Final')
    axs[1, 1].set_xlabel('Posição')
    axs[1, 1].set_ylabel('Probabilidade (%)')
    axs[1, 1].legend()

    # Ajustar layout
    plt.tight_layout()
    plt.show()
    # Adicionar o gráfico ao canvas
    canvas_grafico = FigureCanvasTkAgg(fig, master=frame_grafico)
    canvas_grafico.draw()
    canvas_grafico.get_tk_widget().pack()

def GraficosProtonEletron():
    # Esconder a janela principal
    root.withdraw()

    # Criar uma nova janela para inserir os dados
    janela_probabilidade = tk.Toplevel()
    janela_probabilidade.title("Graficos Proton")
    janela_probabilidade.geometry("600x400")
    root.deiconify()  # Tornar visível a janela principal

    # Campos para inserir os dados
    limite_inf_label = ttk.Label(janela_probabilidade, text="Limite Inferior (nm):")
    limite_inf_label.pack()
    limite_inf_entry = ttk.Entry(janela_probabilidade)
    limite_inf_entry.pack()

    limite_sup_label = ttk.Label(janela_probabilidade, text="Limite Superior (nm):")
    limite_sup_label.pack()
    limite_sup_entry = ttk.Entry(janela_probabilidade)
    limite_sup_entry.pack()

    n1_label = ttk.Label(janela_probabilidade, text="nível de energia inicial:")
    n1_label.pack()
    n1_entry = ttk.Entry(janela_probabilidade)
    n1_entry.pack()

    n2_label = ttk.Label(janela_probabilidade, text="nível de energia final:")
    n2_label.pack()
    n2_entry = ttk.Entry(janela_probabilidade)
    n2_entry.pack()

    l_label = ttk.Label(janela_probabilidade, text="Largura do Poço (nm):")
    l_label.pack()
    l_entry = ttk.Entry(janela_probabilidade)
    l_entry.pack()

    # Frame para o gráfico
    global frame_grafico
    frame_grafico = ttk.Frame(janela_probabilidade)
    frame_grafico.pack(padx=10, pady=10)

    # Botão para plotar os gráficos
    plotar_botao = ttk.Button(janela_probabilidade, text="Plotar Gráficos", command=lambda: plot_grafico(
    limite_inf_entry.get(), limite_sup_entry.get(), n1_entry.get(), n2_entry.get(), l_entry.get()))
    plotar_botao.pack()


def iniciar_simulacao():
    L = float(entry.get())  # Captura o valor de L da entrada
    n = int(mode_entry.get())  # Captura o modo de vibração escolhido

    A = 1  # amplitude da oscilação
    omega = 2 * np.pi / 3  # frequência angular (maior velocidade)
    fps = 30  # quadros por segundo
    t_max = 5  # tempo total da animação

    fig, ax = plt.subplots()
    ax.set_xlim(0, L)
    ax.set_ylim(-1.5 * A, 1.5 * A)
    ax.set_xlabel(f'Posição na Corda (Nível {n})')
    ax.set_title(f'n = {n}')
    line, = ax.plot([], [], lw=2)

    def init():
        line.set_data([], [])
        return line,

    def update(t):
        x = np.linspace(0, L, 1000)
        y = A * np.sin(omega * t) * np.sin(n * np.pi * x / L)
        line.set_data(x, y)
        return line,

    ani = FuncAnimation(fig, update, frames=np.linspace(0, t_max, int(t_max * fps)), init_func=init, blit=True)

    canvas = FigureCanvasTkAgg(fig, master=janela_simulacao)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

def configurar_simulacao():
    global janela_simulacao, entry, mode_entry  # Declara como global

    # Cria a janela de simulação
    janela_simulacao = tk.Toplevel()
    janela_simulacao.title("Simulação de Ondas Animada")
    janela_simulacao.geometry("800x600")

    # Campo de entrada para o comprimento da corda
    entry_label = tk.Label(janela_simulacao, text="Digite o comprimento da corda (L):")
    entry_label.pack()
    entry = tk.Entry(janela_simulacao)
    entry.pack()

    # Campo de entrada para o modo de vibração
    mode_label = tk.Label(janela_simulacao, text="Digite o n:")
    mode_label.pack()
    mode_entry = tk.Entry(janela_simulacao)
    mode_entry.pack()

    # Botão para iniciar a simulação
    iniciar_botao = ttk.Button(janela_simulacao, text="Iniciar Simulação Animada", command=iniciar_simulacao)
    iniciar_botao.pack()
def AbrirOutraJanelaProton():
    # Esconder a janela principal
    root.withdraw()

    # Criar uma nova janela para inserir os dados
    janela_probabilidade = tk.Toplevel()
    janela_probabilidade.title("Próton")
    janela_probabilidade.geometry("1200x600")
    root.deiconify()  # Tornar visível a janela principalroot.deiconify()  # Tornar visível a janela principal
    # Campos para inserir os dados

    limite_inf_label = tk.Label(janela_probabilidade, text="Limite Inferior (nm):")
    limite_inf_label.pack()
    limite_inf_entry = tk.Entry(janela_probabilidade)
    limite_inf_entry.pack()

    limite_sup_label = tk.Label(janela_probabilidade, text="Limite Superior (nm):")
    limite_sup_label.pack()
    limite_sup_entry = tk.Entry(janela_probabilidade)
    limite_sup_entry.pack()

    n1_label = tk.Label(janela_probabilidade, text="nível de energia inicial:")
    n1_label.pack()
    n1_entry = tk.Entry(janela_probabilidade)
    n1_entry.pack()

    n2_label = tk.Label(janela_probabilidade, text="nível de energia final:")
    n2_label.pack()
    n2_entry = tk.Entry(janela_probabilidade)
    n2_entry.pack()

    l_label = tk.Label(janela_probabilidade, text="Largura do Poço (nm):")
    l_label.pack()
    l_entry = tk.Entry(janela_probabilidade)
    l_entry.pack()

# Botão para calcular a probabilidade
    calcular_botao = tk.Button(janela_probabilidade, text="Calcular", command=lambda: calculosProton(
        resultado_text, limite_inf_entry.get(), limite_sup_entry.get(), n1_entry.get(), n2_entry.get(), l_entry.get()))
    calcular_botao.pack()

    # Texto para mostrar o resultado
    resultado_label = tk.Label(janela_probabilidade, text="Resultado:")
    resultado_label.pack()
    resultado_text = tk.Text(janela_probabilidade, height=10, width=50)
    resultado_text.pack()

    
    
    graficos_botao = tk.Button(janela_probabilidade, text="Mostrar Gráficos", command = GraficosProtonEletron)
    graficos_botao.place(x = 1075, y = 150)

    # Botão que abre a janela de configuração da simulação
    configurar_botao = ttk.Button(janela_probabilidade, text="Mostrar Simulação", command=configurar_simulacao)
    configurar_botao.place(x=0, y=150)

    

def AbrirOutraJanelaEletron():
    # Esconder a janela principal
    root.withdraw()

    # Criar uma nova janela para inserir os dados
    janela_probabilidade = tk.Toplevel()
    janela_probabilidade.title("Elétron")
    janela_probabilidade.geometry("1200x800")
    root.deiconify()  # Tornar visível a janela principal

    # Campos para inserir os dados
    limite_inf_label = tk.Label(janela_probabilidade, text="Limite Inferior (nm):")
    limite_inf_label.pack()
    limite_inf_entry = tk.Entry(janela_probabilidade)
    limite_inf_entry.pack()

    limite_sup_label = tk.Label(janela_probabilidade, text="Limite Superior (nm):")
    limite_sup_label.pack()
    limite_sup_entry = tk.Entry(janela_probabilidade)
    limite_sup_entry.pack()

    n1_label = tk.Label(janela_probabilidade, text="nível de energia inicial:")
    n1_label.pack()
    n1_entry = tk.Entry(janela_probabilidade)
    n1_entry.pack()

    n2_label = tk.Label(janela_probabilidade, text="nível de energia final:")
    n2_label.pack()
    n2_entry = tk.Entry(janela_probabilidade)
    n2_entry.pack()

    l_label = tk.Label(janela_probabilidade, text="Largura do Poço (nm):")
    l_label.pack()
    l_entry = tk.Entry(janela_probabilidade)
    l_entry.pack()

    
    
    # Botão para calcular a probabilidade
    calcular_botao = tk.Button(janela_probabilidade, text="Calcular", command=lambda: calculosEletron(
        resultado_text, limite_inf_entry.get(), limite_sup_entry.get(), n1_entry.get(), n2_entry.get(), l_entry.get()))
    calcular_botao.pack()


    resultado_label = tk.Label(janela_probabilidade, text="Resultado:")
    resultado_label.pack()
    resultado_text = tk.Text(janela_probabilidade, height=10, width=50)
    resultado_text.pack()
    
    graficos_botao = tk.Button(janela_probabilidade, text="Mostrar Gráficos", command = GraficosProtonEletron)
    graficos_botao.place(x = 1075, y = 150)

    configurar_botao = ttk.Button(janela_probabilidade, text="Mostrar Simulação", command=configurar_simulacao)
    configurar_botao.place(x=0, y=150)

def botoes():
    #Eletron
    botao = tk.Button(root, text="Confinando um elétron",command=AbrirOutraJanelaEletron,width=20,height=5,bg = "purple", font=("Arial",12))
    botao.place(x=350, y= 200)
    #Próton
    botao2 = tk.Button(root, text="Confinando um próton",command=AbrirOutraJanelaProton,width=20,height=5,bg = "purple", font=("Arial",12))
    botao2.place(x=700, y= 200)
    #Poço de potencial unidimensional
    botao3 = tk.Button(root, text="Poço de potencial unidimensional",command=mostrar_explicacao,width=30,height=5,bg = "purple", font=("Arial",12))
    botao3.place(x=475, y= 50)

# Função para retornar à janela principal
def retornar_janela_principal():
    root.deiconify()  # Tornar visível a janela principal
def PocoDePotencialUnidimensinal():
    botao3 = tk.Button(root, text="Poço de potencial unidimensional",command=mostrar_explicacao,width=30,height=5,bg = "purple", font=("Arial",12))
    botao3.place(x=475, y= 50)
# Função para mostrar uma mensagem explicativa
#def mostrar_explicacao():
    #messagebox.showinfo("Explicação", "O potencial infinito é um cenário físico em que uma partícula está presa por uma barreira de potencial de valor infinito, impossibilitando sua saída. O poço de potencial infinito exemplifica esse conceito na mecânica quântica, onde a função de onda da partícula é nula fora do poço e assume uma forma específica dentro dele, determinada pelas condições de contorno. As energias permitidas para a partícula dentro do poço são quantizadas e relacionadas ao tamanho do poço e ao estado quantizado da partícula. Os autoestados de energia representam os níveis quantizados nos quais a partícula pode existir, sendo soluções da equação de Schrödinger para o poço de potencial infinito. Este modelo é essencial para compreender a quantização de energia, os estados estacionários e as propriedades das funções de onda em sistemas quânticos confinados.")
def mostrar_explicacao():
    # Cria uma nova janela top-level
    janela_explicacao = tk.Toplevel(root)
    janela_explicacao.title("Explicação do Poço de Potencial Unidimensional")
    janela_explicacao.geometry("600x400")

    # Adiciona um texto explicativo
    texto_explicacao = """
    Alunos: Matheus Dourado e João Pedro Sabino

    O potencial infinito é um cenário físico em que uma partícula está presa por uma barreira de potencial de valor infinito, impossibilitando sua saída. O poço de potencial infinito exemplifica esse conceito na mecânica quântica, onde a função de onda da partícula é nula fora do poço e assume uma forma específica dentro dele, determinada pelas condições de contorno. As energias permitidas para a partícula dentro do poço são quantizadas e relacionadas ao tamanho do poço e ao estado quantizado da partícula. Os autoestados de energia representam os níveis quantizados nos quais a partícula pode existir, sendo soluções da equação de Schrödinger para o poço de potencial infinito. Este modelo é essencial para compreender a quantização de energia, os estados estacionários e as propriedades das funções de onda em sistemas quânticos confinados.
    """
    label_explicacao = tk.Label(janela_explicacao, text=texto_explicacao, justify=tk.LEFT, wraplength=550,font=("Times New Roman", 12, "italic"))
    label_explicacao.pack(padx=20, pady=20)
# Criar a janela principal
root = tk.Tk()
root.title("Projeto Simulação Computacional")
root.geometry("1200x600")
root.configure(background="white")
botoes()
# Botão para abrir a janela de cálculo de probabilidade
##abrir_janela_botao.pack()

# Inicia o loop principal da aplicação
root.mainloop()

def probabilidade():

    integral = 0
    precisaoIntegral = 1000
    limiteInf = float(input("Digite o limite inferior do intervalo em nm: "))
    limiteSup = float(input("Digite o limite superior do intervalo em nm: "))
    n = int(input("Digite o nível de energia: "))
    L = float(input("Digite a largura do poço em nm: "))
      
    dx = (limiteInf - limiteSup) / precisaoIntegral  
    
    for i in range(precisaoIntegral):
        x = limiteInf + i * dx
        integral += np.abs(densidade_pontoX(x, n, L))**2 * dx

    porcentagem = integral * 100
    #print(f"sua probabilidade é {abs(porcentagem)} %")
    print(f"A probabilidade é", "{:.2e}".format(abs(porcentagem)) , "%")





constante_normalizacao  = float(input("Digite a constante de normalização: "))
frequencia = float(input("Digite do lado do sen: "))

# Calculando a largura L do poço de potencial
L = 2 / (constante_normalizacao ** 2)  # Convertendo para metros
taxa = float(input("Digite a taxa: "))
x = taxa * L # x = 0.695 * L


probabilidade = (constante_normalizacao ** 2) * np.sin(frequencia * x) ** 2


print("A probabilidade de encontrar o elétron na posição x =", x, "é:", probabilidade)


