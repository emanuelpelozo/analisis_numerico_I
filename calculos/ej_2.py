import numpy as np
import matplotlib.pyplot as plt
from ej_1 import pasos_h, m, k, y0,v0, INTERVALO_T

nuevo_lamda = 750
lambda_ = nuevo_lamda
def valor_de_c(tiempo):
    """
    Devuelve el valor de c y c' en funcion del tiempo, en forma de una tupla (c , c').
    """

    if (tiempo <= 1.1 and tiempo >1.0 ) :
        return (tiempo-1,1)
    
    elif (tiempo<1.4 and tiempo>=1.3):
        return (tiempo-1.3,-1)
    
    elif (tiempo<1.3 and tiempo>1.1):
        return(0.1,0)
    
    else:
        return (0,0)
    

def calculo_q_c_dinamico(u,v,paso,punto_tiempo):
    """
    Calculamos las q's
    """
    c, c_prima = valor_de_c(punto_tiempo)
    q1u = paso*v
    q1v = paso*((1/m)*(k*(c-u)+lambda_*(c_prima-v)))

    q2u = paso*(v+q1v)
    q2v = paso*((1/m)*(k*(c-(u+q1u))+lambda_*(c_prima-(v+q1v))))

    return q1u,q2u,q1v,q2v

def runge_kutta_c_dinamico(funcion, u0, v0, paso_h, intervalo_t):
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

        q1u,q2u,q1v,q2v = calculo_q_c_dinamico(u,v,paso_h, t)

        u += 0.5*(q1u+q2u)
        v += 0.5*(q1v+q2v)

    return tiempo, valores_u, valores_v

# Gráfico comparativo
plt.figure(figsize=(10, 6))

for paso_h in pasos_h:
    tiempo, valores_y, _ = runge_kutta_c_dinamico(calculo_q_c_dinamico, y0, v0, paso_h, INTERVALO_T)
    plt.plot(tiempo, valores_y, label=f"Paso h={paso_h}")

plt.title("Método de Runge Kutta - Comparación de Pasos")
plt.xlabel("Tiempo (s)")
plt.ylabel("Posición (m)")
plt.legend()
plt.grid(True)
plt.show()
