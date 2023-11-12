import numpy as np
import matplotlib.pyplot as plt

PADRON = 99444
# Variables a usar:
m = PADRON / 200  # kg
k = 25000  # N/m
lambda_ = 0  # Ns/m
c = 0.1  # m
c_prima = 0 # estoy asumiendo c' como 0
paso_h = 0.005  # segundos
pasos_h = [0.005, 0.001, 0.01]  # Lista de diferentes tamaños de paso en segundos
INTERVALO_T = 5  # segundos

# Aplico cambio de variable:
def sistema_cambio_de_variable(u, v):
    u_prima = v
    v_prima = (k / m) * (c - u) + (lambda_ / m) * (c_prima- v) 
    return u_prima, v_prima

# Euler explicito
def euler_explicito(funcion, u0, v0, paso_h, intervalo_t):
    # Creamos vector tiempo desde 0 hasta intervalo_t, con un paso_h. 
    # Aca guardamos los puntos en el tiempo donde vamos a evaluar la solución.
    tiempo = np.arange(0, intervalo_t, paso_h)

    valores_u = []
    valores_v = []

    # seteamos condiciones iniciales
    u = u0
    v = v0

    # Recorremos los puntos en el tiempo.
    for t in tiempo:
        # Guardo el valor obtenido para el paso actual en los vectores.
        valores_u.append(u)
        valores_v.append(v)

        u_prima, v_prima = funcion(u, v)

        u += u_prima * paso_h
        v += v_prima * paso_h

    return tiempo, valores_u, valores_v


## Trato de resolver el problema:

y0 = 0  # Condición inicial para y
v0 = 0  # Condición inicial para v


# tiempo, valores_u, valores_v = euler_explicito(sistema_cambio_de_variable, y0, v0, paso_h, INTERVALO_T)

# # Graficar resultados
# plt.plot(tiempo, valores_u, label='Posición (u)')
# #plt.plot(tiempo, valores_v, label='Velocidad (v)')
# plt.xlabel('Tiempo (s)')
# plt.ylabel('Magnitud')
# plt.legend()
# plt.show()


# Gráfico comparativo
plt.figure(figsize=(10, 6))

for paso_h in pasos_h:
    tiempo, valores_y, _ = euler_explicito(sistema_cambio_de_variable, y0, v0, paso_h, INTERVALO_T)
    plt.plot(tiempo, valores_y, label=f"Paso h={paso_h}")

plt.title("Método de Euler Explícito - Comparación de Pasos")
plt.xlabel("Tiempo (s)")
plt.ylabel("Posición (m)")
plt.legend()
plt.grid(True)
plt.show()
