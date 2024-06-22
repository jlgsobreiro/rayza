from math import exp, log
import matplotlib.pyplot as plt


def IDCPH(T0, T, deltaA, deltaB, deltaC, deltaD):
    return deltaA * (T - T0) + deltaB * log(T / T0) + deltaC + deltaD * (T - T0) ** 2


def IDCPS(T0, T, deltaA, deltaB, deltaC, deltaD):
    return deltaA * (T - T0) + deltaB * (T / T0) + deltaC + deltaD * (T - T0)


print('''13.48. O craqueamento de propano é uma rota para a produção de olefinas leves. Supona que duas reaçoes de
craqueamento ocorram em um reator contínuo em regime estacionário:
C3H8(g) -> C3H6(g) + H2(g)      (I)
C3H8(g) -> C2H4(g) + CH4(g)     (II)
Calcule a composição do produto se as duas reações alcançam o equilibrio a 1,2 bar e
(a) 750K; (b) 1OOOK; (c) 1250K''')

'''13.48
C3H8(g) -> C3H6(g) + H2(g)  (I)
C3H8(g) -> C2H4(g) + CH4(g) (II)
TO = 298.15K
PO = 1 (bar)
T = 750K
P = 1.2(bar)
'''
# # 1 = C3H8(g)
#
# deltaH0f1 := -104680  #J/mol
# deltaG0f1 := -24290   #J/mol
#
# # 2 = C3H6(g)
#
# deltaH0f2 := -19710  #J/mol
# deltaG0f2 := -62205  #J/mol
#
# # 3 = H2(g)
#
# deltaH0f3 := 0  #J/mol
# deltaG0f3 := 0  #J/mol
#
# # 4 = C2H4 (g)
#
# deltaH0f4 := 52510  #J/mol
# deltaG0f4 := 68460   #J/mol
#
# # 5 = CH4 (g)
#
# deltaH0f5 := -74520  #J/mol
# deltaG0f5 := -50460  #J/mol

R = 8.31  #J/mol.K

T0 = 298.15  #K
T_list = range(750, 1300, 50)  #K
P0 = 1  #(bar)
P = 1.2  #(bar)

dict_elements = {
    'C3H8': {"deltah0": -104680, "deltag0": -24290},
    'C3H6': {"deltah0": 19710, "deltag0": 62205},
    'H2': {"deltah0": 0, "deltag0": 0},
    'C2H4': {"deltah0": 52510, "deltag0": 68460},
    'CH4': {"deltah0": -74520, "deltag0": -50460}
}

tabela_C1 = {
    "C3H8": {"Vi": -1, "Ai": 1.213, "Bi": 28.785 * 10 ** -3, "Ci": -8.824 * 10 ** -6, "Di": 0},
    "C3H6": {"Vi": 1, "Ai": 1.637, "Bi": 22.706 * 10 ** -3, "Ci": -6.915 * 10 ** -6, "Di": 0},
    "H2": {"Vi": 1, "Ai": 3.249, "Bi": 0.422 * 10 ** -3, "Ci": 0, "Di": 0.083 * 10 ** 5},
    "C2H4": {"Vi": 1, "Ai": 1.424, "Bi": 14.394 * 10 ** -3, "Ci": -4.392 * 10 ** -6, "Di": 0},
    "CH4": {"Vi": 1, "Ai": 1.702, "Bi": 9.081 * 10 ** -3, "Ci": -2.164 * 10 ** -6, "Di": 0},
}


def delta_H0I(a, b, c):
    return - dict_elements[a]["deltah0"] + dict_elements[b]["deltah0"] + dict_elements[c]["deltah0"]


def delta_G0I(a, b, c):
    return - dict_elements[a]["deltag0"] + dict_elements[b]["deltag0"] + dict_elements[c]["deltag0"]


def deltaAI(a, b, c):
    return - tabela_C1[a]["Ai"] + tabela_C1[b]["Ai"] + tabela_C1[c]["Ai"]


def deltaBI(a, b, c):
    return - tabela_C1[a]["Bi"] + tabela_C1[b]["Bi"] + tabela_C1[c]["Bi"]


def deltaCI(a, b, c):
    return - tabela_C1[a]["Ci"] + tabela_C1[b]["Ci"] + tabela_C1[c]["Ci"]


def deltaDI(a, b, c):
    return - tabela_C1[a]["Di"] + tabela_C1[b]["Di"] + tabela_C1[c]["Di"]


def tau(T0, T):
    return T / T0


def KI0(deltaG0I, T0):
    return exp(-deltaG0I / (R * T0))


def KI1(deltaH0I, T0, T):
    return exp((deltaH0I / (R * T0)) * (1 - T0 / T))


def KI2(T0, T, dAI, dBI, dCI, dDI):
    tau0 = tau(T0, T)
    result = exp(
        dAI * (log(tau0) - ((tau0 - 1) / tau0)) +
        (1 / 2) * dBI * T0 * (((tau0 - 1) ** 2) / tau0) +
        (1 / 6) * (dCI * T0 ** 2) * ((((tau0 - 1) ** 2) * (tau0 + 2)) / tau0) +
        (1 / 2) * (dDI / T0 ** 2) * (((tau0 - 1) ** 2) / tau0 ** 2))
    return result


def K(k0, k1, k2):
    return k0 * k1 * k2


def V(a, b, c):
    return tabela_C1[a]["Vi"] + tabela_C1[b]["Vi"] + tabela_C1[c]["Vi"]


def n0(a, b, c):
    return - tabela_C1[a]["Vi"] + tabela_C1[b]["Vi"] + tabela_C1[c]["Vi"]


reactions = [("C3H8", "C3H6", "H2"), ("C3H8", "C2H4", "CH4")]
Ks = []
for reaction in reactions:
    for T in T_list:
        print(f'\n{reaction[0]} -> {reaction[1]} + {reaction[2]}')
        print(f'T = {T}K')
        dH0I = delta_H0I(reaction[0], reaction[1], reaction[2])
        dG0I = delta_G0I(reaction[0], reaction[1], reaction[2])
        print(f'deltaH0I = {dH0I} J/mol')
        print(f'deltaG0I = {dG0I} J/mol')
        dAI = deltaAI(reaction[0], reaction[1], reaction[2])
        dBI = deltaBI(reaction[0], reaction[1], reaction[2])
        dCI = deltaCI(reaction[0], reaction[1], reaction[2])
        dDI = deltaDI(reaction[0], reaction[1], reaction[2])
        print(f'deltaAI = {dAI}')
        print(f'deltaBI = {dBI}')
        print(f'deltaCI = {dCI}')
        print(f'deltaDI = {dDI}')
        tau0 = tau(T0, T)
        print(f'tau0 = {tau0}')
        K0 = KI0(dG0I, T0)
        print(f'KI0 = {K0}')
        K1 = KI1(dH0I, T0, T)
        print(f'KI1 = {K1}')
        K2 = KI2(T0, T, dAI, dBI, dCI, dDI)
        print(f'KI2 = {K2}')
        KI = K(K0, K1, K2)
        Ks.append(KI)
        print(f'KI = {KI}')
        v = V(reaction[0], reaction[1], reaction[2])
        print(f'v = {v}')
        n = n0(reaction[0], reaction[1], reaction[2])
        print(f'n0 = {n}')

KI_list = Ks[:len(T_list)]
KII_list = Ks[len(T_list):]
epsilonI = 0.5
epsilonII = 0.5
p = P0 / P

y_list = []

for i in range(len(T_list)):
    print(f'\nREAÇÃO {i + 1}')

    KI = KI_list[i//2]
    KII = KII_list[i//2]
    print(f'\nT = {T_list[i]}K')
    print(f'KI = {KI}')
    print(f'KII = {KII}')


    def equation(k):
        eq = 1/((1 / (1 + epsilon_I + epsilon_II)) ** 2 - (k * p * (1 - epsilon_I - epsilon_II) / (1 + epsilon_I + epsilon_II)))
        return eq

    def equation_derivative(k):
        eq = (epsilon_I**2 / (1 + epsilon_I + epsilon_II) ** 2) - ((k * p * (1 - epsilon_I - epsilon_II)) / (1 + epsilon_I + epsilon_II))
        return eq

    def eq2(k):
        return ((k*p*((1-epsilonI-epsilonII)/(1+epsilonI+epsilonII)))**(1/2)) * (1+epsilonI+epsilonII)

    def newton_raphson(f, fdx, tol, max_iter=1):
        x = None
        for i in range(max_iter):
            x = tol - f(tol) / fdx(tol)
        return x

    # epsilon_I = newton_raphson(equation, equation_derivative, KI)
    # epsilon_II = newton_raphson(equation, equation_derivative, KII)
    epsilon_I = eq2(KI)
    epsilon_II = eq2(KII)
    y1 = (1 - epsilon_I - epsilon_II) / (1 + epsilon_I + epsilon_II)
    y2 = epsilon_I / (1 + epsilon_I + epsilon_II)
    y3 = epsilon_I / (1 + epsilon_I + epsilon_II)
    y4 = epsilon_II / (1 + epsilon_I + epsilon_II)
    y5 = epsilon_II / (1 + epsilon_I + epsilon_II)
    print(f"epsilon_I = {epsilon_I}")
    print(f"epsilon_II = {epsilon_II}")
    print(f"y1 = {y1}")
    print(f"y2 = {y2}")
    print(f"y3 = {y3}")
    print(f"y4 = {y4}")
    print(f"y5 = {y5}")
    y_list.append([y1, y2, y3, y4, y5])

print(f'\n\nRESULTADO FINAL')
for i in range(len(y_list)):
    print(f'\nT = {T_list[i]}K')
    print(f'y1 = {y_list[i][0]}')
    print(f'y2 = {y_list[i][1]}')
    print(f'y3 = {y_list[i][2]}')
    print(f'y4 = {y_list[i][3]}')
    print(f'y5 = {y_list[i][4]}')

# Resultado final
plt.plot(T_list, [y[0] for y in y_list], label='y1')
plt.plot(T_list, [y[1] for y in y_list], label='y2')
plt.plot(T_list, [y[2] for y in y_list], label='y3')
plt.plot(T_list, [y[3] for y in y_list], label='y4')
plt.plot(T_list, [y[4] for y in y_list], label='y5')
plt.legend()
plt.show()
