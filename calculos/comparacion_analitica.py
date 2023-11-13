import numpy as np
import matplotlib.pyplot as plt
import constantes as cts
from ej_1 import euler_explicito, sistema_cambio_de_variable
from scipy.signal import find_peaks

def ejecutar_euler_explicito(y0, v0, paso_h, intervalo):
    return euler_explicito(sistema_cambio_de_variable, y0, v0, paso_h, intervalo)


# Calcular la frecuencia natural de oscilación
frecuencia_calculada = np.sqrt(cts.k / cts.m)

# Resuelvo con euler explicito
tiempo, valores_y, valores_v = ejecutar_euler_explicito(cts.y0, cts.v0, 0.0005, cts.INTERVALO_T)
# Encuentra los índices de los máximos de los puntos obtenidos numericamente
indices_maximos, _ = find_peaks(valores_y)


# Graficar las oscilaciones y marcar los máximos
plt.plot(tiempo, valores_y, label='Posición (y)')


# Marca los máximos en el gráfico
for indice in indices_maximos:
    plt.axvline(tiempo[indice], color='r', linestyle='--', alpha=0.5)

plt.xlabel('Tiempo (s)')
plt.ylabel('Posición (m)')
plt.title('Oscilaciones del sistema y Máximos')
plt.legend()
plt.show()

# Calcula las diferencias de tiempo entre extremos consecutivos
diferencias_tiempo = np.diff(tiempo[indices_maximos])

# Calcula el período promedio
periodo_promedio = np.mean(diferencias_tiempo)

# Calcula la frecuencia angular de forma experimental en base al periodo obtenido
frecuencia_angular_experimental = 2 * np.pi / periodo_promedio

print("Período promedio:", periodo_promedio)
print("Frecuencia angular experimental:", frecuencia_angular_experimental)
print(f"Frecuencia natural de oscilación calculada: {frecuencia_calculada} rad/s")
print(f"Relacion: frecuenca_natural/frecuencia_analitica: {frecuencia_calculada/frecuencia_angular_experimental} ")

