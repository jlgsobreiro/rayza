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


def gauss_jordan(A, b):
    n = len(A)
    M = A
    for i in range(n):
        M[i].append(b[i])
    for i in range(n):
        M[i] = [m_ij / M[i][i] for m_ij in M[i]]
        for j in range(n):
            if i != j:
                M[j] = [M[j][k] - M[i][k]*M[j][i] for k in range(n+1)]
    return [M[i][n] for i in range(n)]


def jacobian(F, x, h=1e-5):
    n = len(x)
    J = [[0]*n for _ in range(n)]
    for i in range(n):
        x[i] += h
        F_plus_h = F(x)
        x[i] -= 2*h
        F_minus_h = F(x)
        x[i] += h
        J[i] = [(F_plus_h_j - F_minus_h_j) / (2*h) for F_plus_h_j, F_minus_h_j in zip(F_plus_h, F_minus_h)]
    return J


def newton_raphson(F, x, eps=5e-5):
    F_value = F(x)
    F_norm = sum(F_value_i ** 2 for F_value_i in F_value) ** 0.5  # l2 norm of vector
    iteration_counter = 0
    while abs(F_norm) > eps and iteration_counter < 100:
        J = jacobian(F, x)
        delta = gauss_jordan(J, [-F_value_i for F_value_i in F_value])
        x = [x_i + delta_i for x_i, delta_i in zip(x, delta)]
        F_value = F(x)
        F_norm = sum(F_value_i ** 2 for F_value_i in F_value) ** 0.5
        iteration_counter += 1

    if abs(F_norm) > eps:
        iteration_counter = -1
    return x, iteration_counter


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
epsilon_I = 0.5
epsilon_II = 0.5
p = P0 / P

y_list = []

for i in range(len(T_list)):
    print(f'\nREAÇÃO {i + 1}')

    KI = KI_list[i//2]
    KII = KII_list[i//2]
    print(f'\nT = {T_list[i]}K')
    print(f'KI = {KI}')
    print(f'KII = {KII}')


    def equations(vars):
        epsilon_I, epsilon_II = vars
        eq1 = (epsilon_I / (1 + epsilon_I + epsilon_II)) ** 2 - (
                    KI * p * (1 - epsilon_I - epsilon_II) / (1 + epsilon_I + epsilon_II))
        eq2 = (epsilon_II / (1 + epsilon_I + epsilon_II)) ** 2 - (
                    KII * p * (1 - epsilon_I - epsilon_II) / (1 + epsilon_I + epsilon_II))
        return [eq1, eq2]


    solution, no_iterations = newton_raphson(equations, [epsilon_I, epsilon_II])
    epsilon_I, epsilon_II = solution
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
