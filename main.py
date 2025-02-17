import numpy as np
import matplotlib.pyplot as plt

a = 0.1382
b = 3.19 * 10 ** (-5)
R = 8.314
Pnas = 3664186.998


def P(V, T, R=R, a=a, b=b):
    return R * T / (V - b) - a / V ** 2


temp = np.array([-140, -130, -120, -110, -100])
Temp = temp + 273.15
colors = np.array(['green', 'yellow', 'red', 'blue', 'orange'])
for i in range(5):
    T = Temp[i]
    color = colors[i]
    V = np.linspace(b + 10 ** (-5), 3 * 10 ** (-4), 2000)
    P_values = P(V, T)
    plt.plot(V, P_values, label=f"T = {temp[i]}°C", color=color)
plt.xlabel("V, м ** 3 / моль")
plt.ylabel("P, Па")
plt.legend()
plt.show()
# Начиная с температуры -120°C температура превышает критическую
T = -130 + 273.15


def find_min(a, b, f, eps):
    while b - a > 2 * eps:
        mid1 = a + (b - a) / 3
        mid2 = b - (b - a) / 3
        if f(mid1) < f(mid2):
            b = mid2
        else:
            a = mid1
    return (a + b) / 2


def find_max(a, b, f, eps):
    while b - a > 2 * eps:
        mid1 = a + (b - a) / 3
        mid2 = b - (b - a) / 3
        if f(mid1) > f(mid2):
            b = mid2
        else:
            a = mid1
    return (a + b) / 2


def calc_length(a, b, f, n=10000):
    h = (b - a) / n
    x = a
    y = f(x)
    length = 0
    while x < b:
        x += h
        y_next = f(x)
        length += np.sqrt(h ** 2 + (y_next - y) ** 2)
        y = y_next
    return length


def solve_equation(a, b, f, eps):
    while b - a >= 2 * eps:
        c = (a + b) / 2
        if f(c) * f(a) < 0:
            b = c
        else:
            a = c
    return c


def square(a, b, f, n=10000):
    h = (b - a) / n
    x = a
    y = f(x)
    area = 0
    while x < b:
        x += h
        y_next = f(x)
        area += 0.5 * (y + y_next) * h
        y = y_next
    return area


V_1 = find_min(b + 10 ** (-5), 10 ** (-4), lambda x: P(x, T=T), 10 ** (-10))
V_2 = find_max(10 ** (-4), 3 * 10 ** (-4), lambda x: P(x, T=T), 10 ** (-10))
print(f'Локальный минимум: {V_1}')
print(f'Локальный максимум: {V_2}')
print(f'Длина кривой запрещённой зоны: {calc_length(V_1, V_2, lambda x: P(x, T=T))}')

V1 = solve_equation(5 * 10 ** (-5), 9 * 10 ** (-5), lambda x: P(x, T) - Pnas, 10 ** (-10))
V_c = solve_equation(9 * 10 ** (-5), 11 * 10 ** (-5), lambda x: P(x, T) - Pnas, 10 ** (-10))
V_g = solve_equation(18 * 10 ** (-5), 21 * 10 ** (-5), lambda x: P(x, T) - Pnas, 10 ** (-10))
print(f'Корни уравнения:\n{V1}\n{V_c}\n{V_g}')
print(f'Первый интеграл: {square(V1, V_g, lambda x: P(x, T=T))}')
print(f'Второй интеграл: {square(V1, V_g, lambda x: Pnas)}')
