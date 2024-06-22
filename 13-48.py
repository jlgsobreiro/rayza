import matplotlib.pyplot as plt

# Função exp e log manualmente
def exp(x, tol=1e-10):
    term = 1
    sum = 1
    n = 1
    while abs(term) > tol:
        term *= x / n
        sum += term
        n += 1
    return sum

def log(x, tol=1e-10):
    if x <= 0:
        raise ValueError("Logarithm only defined for positive values")
    n = 0
    while x > 2:
        x /= 2
        n += 1
    x -= 1
    term = x
    sum = x
    i = 2
    while abs(term) > tol:
        term *= -x * (i - 1) / i
        sum += term
        i += 1
    return sum + n * 0.69314718056

# Constantes e dados
R = 8.31  # J/mol.K
T0 = 298.15  # K
T_list = range(750, 1300, 50)  # K
P0 = 1  # bar
P = 1.2  # bar

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


    # Método de Newton-Raphson
    def jacobian(epsilons, KI, KII, p):
        epsilon_I, epsilon_II = epsilons
        common_denominator = (1 + epsilon_I + epsilon_II) ** 3
        d_eq1_d_epsilon_I = (2 * epsilon_I * (1 + epsilon_I + epsilon_II) - epsilon_I ** 2 - epsilon_II ** 2 - 1) / common_denominator - (KI * p * (-1)) / (1 + epsilon_I + epsilon_II) ** 2
        d_eq1_d_epsilon_II = (2 * epsilon_I * (1 + epsilon_I + epsilon_II) - epsilon_I ** 2 - epsilon_II ** 2 - 1) / common_denominator
        d_eq2_d_epsilon_I = (2 * epsilon_II * (1 + epsilon_I + epsilon_II) - epsilon_I ** 2 - epsilon_II ** 2 - 1) / common_denominator
        d_eq2_d_epsilon_II = (2 * epsilon_II * (1 + epsilon_I + epsilon_II) - epsilon_I ** 2 - epsilon_II ** 2 - 1) / common_denominator - (KII * p * (-1)) / (1 + epsilon_I + epsilon_II) ** 2
        return [[d_eq1_d_epsilon_I, d_eq1_d_epsilon_II],
                [d_eq2_d_epsilon_I, d_eq2_d_epsilon_II]]


    def newton_raphson(equations, initial_guess, KI, KII, p, tol=1e-6, max_iterations=3):
        epsilon_I, epsilon_II = initial_guess
        for _ in range(max_iterations):
            F = equations([epsilon_I, epsilon_II])
            J = jacobian([epsilon_I, epsilon_II], KI, KII, p)
            # Resolução manual do sistema linear
            det_J = J[0][0] * J[1][1] - J[0][1] * J[1][0]
            if det_J == 0:
                raise ValueError("Jacobian determinant is zero")
            delta_I = (F[0] * J[1][1] - F[1] * J[0][1]) / det_J
            delta_II = (J[0][0] * F[1] - J[1][0] * F[0]) / det_J
            epsilon_I, epsilon_II = epsilon_I - delta_I, epsilon_II - delta_II
            if abs(delta_I) < tol and abs(delta_II) < tol:
                break
        return epsilon_I, epsilon_II


    epsilon_I, epsilon_II = newton_raphson(equations, [epsilon_I, epsilon_II], KI, KII, p)
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
