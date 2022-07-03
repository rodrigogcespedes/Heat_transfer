import math
import numpy as np
import control as co
import sympy
import PySimpleGUI as sg
from matplotlib import pyplot as plt

# Const
s = sympy.symbols('s')
rec = np.linspace(0, 0.001, 1000)
rec2 = np.linspace(0, 0.00005, 1000)
ap1 = 0
es = 0
mp1 = 0
cp1 = 0
mf = 0
cf = 0
k = 0
qE = 0
sistema = "f" # puede ser "f" para analizar la función de transferencia del fluido o "p1" para analizar la función de transferencia de la pared


def gui():
    layout_data = [
        [sg.Text("Ingrese los datos del experimento")],
        [sg.Text('Flujo de enegía entregado por la resistencia:')],
        [sg.Input(key='-qE-')],
        [sg.Text('Area de la pared:')],
        [sg.Input(key='-ap1-')],
        [sg.Text('Espesor de la pared:')],
        [sg.Input(key='-es-')],
        [sg.Text('Coeficiente de conductividad de la pared:')],
        [sg.Input(key='-k-')],
        [sg.Text('Masa de la pared:')],
        [sg.Input(key='-mp1-')],
        [sg.Text('Calor específico de la pared:')],
        [sg.Input(key='-cp1-')],
        [sg.Text('Masa del fluido:')],
        [sg.Input(key='-mf-')],
        [sg.Text('Calor específico del fluido:')],
        [sg.Input(key='-cf-')],
        [sg.Button('Aceptar')]
    ]

    window_data = sg.Window('Sistema Termico', layout_data)

    while True:
        event, values = window_data.read()
        if event == sg.WINDOW_CLOSED:
            exit()
        if event == 'Aceptar':
            ap1 = float(values['-ap1-'])
            es = float(values['-es-'])
            mp1 = float(values['-mp1-'])
            cp1 = float(values['-cp1-'])
            mf = float(values['-mf-'])
            cf = float(values['-cf-'])
            k = float(values['-k-'])
            qE = float(values['-qE-'])
            break

    window_data.close()

    # Primitivos Ej1
    # ap1 = 2
    # es = 0.2
    # mp1 = 800
    # cp1 = 900
    # mf = 9.6328
    # cf = 1005
    # k = 0.75
    # qE = 2000

    # Primitivos Ej2
    # ap1 = 1
    # es = 0.01
    # mp1 = 25
    # cp1 = 840
    # mf = 941.2
    # cf = 4180
    # k = 1.05
    # qE = 2000

    # Procesados
    Rp1 = es / k
    Cf = cf * mf
    Cp1 = cp1 * mp1

    # Funcionp1
    H1 = -1 * ap1 * (Cf + Cp1) / (Rp1 * Cf * Cp1)
    Tp1 = co.tf([qE * ap1 / (Rp1 * Cp1 * Cf)], [1, -H1, 0])
    ss = co.tf2ss(Tp1)
    FCCp1, mat1 = co.canonical_form(ss, form='reachable')
    FCOp1, mat2 = co.canonical_form(ss, form='observable')
    zerosp1 = co.zero(Tp1)
    polosp1 = co.pole(Tp1)

    outEU = qE * ap1 / (Rp1 * Cp1 * Cf * (rec ** 3 - H1 * rec**2))
    kpp1 = sympy.limit(qE * ap1 / (Rp1 * Cp1 * Cf * (s ** 2 - H1 * s)), s, 0)
    errkpp1 = 1 / (1 + kpp1)

    outRU = qE * ap1 / (Rp1 * Cp1 * Cf * (rec ** 4 - H1 * rec ** 3))
    kvp1 = sympy.limit(qE * ap1 * s / (Rp1 * Cp1 * Cf * (s ** 2 - H1 * s)), s, 0)
    errkvp1 = 1 / kvp1

    outPU = qE * ap1 / (Rp1 * Cp1 * Cf * (rec ** 5 - H1 * rec ** 4))
    kap1 = sympy.limit(qE * ap1 * s ** 2 / (Rp1 * Cp1 * Cf * (s ** 2 - H1 * s)), s, 0)
    errkap1 = 1 / kap1

    t1p1 = np.array(range(500))
    x1p1 = np.zeros(len(t1p1))
    for i in range(len(t1p1)):
        x1p1[i] = math.exp(0 * t1p1[i])

    t2p1 = np.array(t1p1)
    x2p1 = np.zeros(len(t2p1))
    for j in range(len(t2p1)):
        x2p1[j] = math.exp(t2p1[j] * H1)

    # Funcionf
    F = ap1 / Rp1
    C = F ** 2 / (Cp1 * Cf)
    B = F * (2 / Cf + 1 / Cp1)
    Tf = co.tf([qE*Cp1/(Cp1*Cf), qE*F/(Cp1*Cf)], [1, B, C])
    ssf = co.tf2ss(Tf)
    FCCf, mat1f = co.canonical_form(ssf, form='reachable')
    FCOf, mat2f = co.canonical_form(ssf, form='observable')
    zerosf = co.zero(Tf)
    polosf = co.pole(Tf)

    outEUf = ((qE*Cp1*rec+F)/(Cp1*Cf))/(rec**3+B*rec**2+C*rec)
    kpf = sympy.limit(((qE*Cp1*s+F)/(Cp1*Cf))/(s**2+B*s+C), s, 0)
    errkpf = 1 / (1 + kpf)

    outRUf = ((qE*Cp1*rec+F)/(Cp1*Cf))/(rec**4+B*rec**3+C*rec**2)
    kvf = sympy.limit(((qE*Cp1*s+F)*s/(Cp1*Cf))/(s**2+B*s+C), s, 0)
    errkvf = 1 / kvf

    outPUf = ((qE * Cp1 * rec + F) / (Cp1 * Cf)) / (rec ** 5 + B * rec ** 4 + C * rec ** 3)
    kaf = sympy.limit(((qE*Cp1*s+F)*s**2/(Cp1*Cf))/(s**2+B*s+C), s, 0)
    errkaf = 1 / kaf

    B = F*(2/Cf+1/Cp1)
    C = F**2 /(Cf*Cp1)
    escalar = qE/(Cf*Cp1)
    nyq1 = escalar * (C*F+rec2**2 * (Cp1*B-F))/(rec2**4+rec2**2*(B-2*C)+C**2)
    nyq2 = escalar*(-1*Cp1*rec2**3+rec2*(Cp1*C-F*B))/(rec2**4+rec2**2*(B-2*C)+C**2)

    plt.subplots_adjust(top=0.949, bottom=0.083, left=0.049, right=0.989, hspace=0.217, wspace=0.291)

    if sistema == "p1":
        # Plot p1
        plt.subplot(241)
        plt.plot(rec, outEU)
        plt.title("Salida Escalon Unitario")

        plt.subplot(242)
        plt.plot(rec, outRU)
        plt.title("Salida Rampa Unitaria")

        plt.subplot(243)
        plt.plot(rec, outPU)
        plt.title("Salida Parabola Unitaria")

        plt.subplot(244)
        plt.axis('off')

        texto = "Funcion de transferencia:\n"\
                + str(Tp1) + "\n"\
                + "Polos de la función: " + str(polosp1).replace(".        ", "") + "\n" \
                + "Ceros de la función: " + str(zerosp1) + "\n" \
                + "Estabilidad: Criticamente estable\n" \
                + "Error Escalon unitario: "\
                + str(errkpp1) + "\n"\
                + "Error Rampa unitaria: " \
                + str(errkvp1) + "\n"\
                + "Error Parabola unitaria: " \
                + str(errkap1) + "\n"\
                + "FCC: \n" \
                + str(FCCp1).replace("\n\n", "\n") + "\n" \
                + "FCO: \n" \
                + str(FCOp1).replace("\n\n", "\n") + "\n" \
                + "FCD: \n" \
                + "A = [[0 0]\n" \
                  "    [0 "+str(H1)+"]]\n" \
                + "B = [[1]\n" \
                  "    [1]]\n" \
                + "C = ["+str(qE/(Cp1+Cf))+" "+str(-qE/(Cp1+Cf))+"]\n"

        plt.text(-0.2, 1.0, texto, fontsize=11, va="top", wrap=True)

        plt.subplot(245)
        plt.plot(t1p1, x1p1, t2p1, x2p1)
        plt.xlabel("tiempo")
        plt.ylabel("x\'1, x\'2")
        plt.legend(['x\'1', 'x\'2'])
        plt.title("Representacion de funciones")

        plt.subplot(246)
        plt.plot(x1p1, x2p1)
        plt.xlabel("x\'1")
        plt.ylabel("x\'2")
        plt.legend(['x\'1', 'x\'2'])
        plt.title("Plano de Fase")

        plt.subplot(247)
        LCp1 = co.nyquist(Tp1 * 1, np.linspace(-10000000, 100000000, 100000))
        plt.plot(LCp1)
        plt.title("Imagen de Nyquist")

    else:
        # Plot f
        plt.subplot(241)
        plt.plot(rec, outEUf)
        plt.title("Salida Escalon Unitario")

        plt.subplot(242)
        plt.plot(rec, outRUf)
        plt.title("Salida Rampa Unitaria")

        plt.subplot(243)
        plt.plot(rec, outPUf)
        plt.title("Salida Parabola Unitaria")

        plt.subplot(244)
        plt.axis('off')

        textof = "Funcion de transferencia:\n"\
                + str(Tf) + "\n"\
                + "Polos de la función: " + str(polosf).replace(".        ", "") + "\n" \
                + "Ceros de la función: " + str(zerosf) + "\n" \
                + "Estabilidad: estable sobreamortiguado\n" \
                + "Error Escalon unitario: "\
                + str(errkpf) + "\n"\
                + "Error Rampa unitaria: " \
                + str(errkvf) + "\n"\
                + "Error Parabola unitaria: " \
                + str(errkaf) + "\n"\
                + "FCC: \n" \
                + str(FCCf).replace("\n\n", "\n") + "\n" \
                + "FCO: \n" \
                + str(FCOf).replace("\n\n", "\n") + "\n" \
                + "FCD: \n" \
                + "A = [["+str((F*(Cf-2*Cp1))/(2*Cf*Cp1) + math.sqrt(((Cf-2*Cp1)**2) /(Cf**2 * Cp1**2)-4/(Cf*Cp1))*F/2)+" 0]\n" \
                + "    [0 " + str((F*(Cf-2*Cp1))/(2*Cf*Cp1) - math.sqrt(((Cf-2*Cp1)**2) /(Cf**2 * Cp1**2)-4/(Cf*Cp1))*F/2) + "]]\n" \
                + "B = [[1]\n" \
                + "    [1]]\n" \
                + "C = [" + str(qE**2 / math.sqrt(((Cf-2*Cp1)**2) / (Cf**2 * Cp1**2)-4/(Cf*Cp1))) + " " + str(-qE**2 / math.sqrt(((Cf-2*Cp1)**2) / (Cf**2 * Cp1**2)-4/(Cf*Cp1))) + "]\n"

        plt.text(-0.2, 1.0, textof, fontsize=11, va="top", wrap=True)

        plt.subplot(246)
        plt.plot(nyq1, nyq2)
        plt.title("Camino de Nyquist")

        plt.subplot(247)
        LCf = co.nyquist(Tf * 1, np.linspace(-10000000, 100000000, 100000))
        plt.plot(LCf)
        plt.title("Imagen de Nyquist")

    plt.show()


if __name__ == '__main__':
    gui()
