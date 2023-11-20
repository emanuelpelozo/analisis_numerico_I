import numpy as np
import matplotlib.pyplot as plt
from math import sqrt as raiz

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



# ----------------------------- EULER EXPLICITO --------------------------------


# Aplico cambio de variable:
def sistema_cambio_de_variable(u, v):
    u_prima = v
    v_prima = (k / m) * (c - u) + (lambda_ / m) * (c_prima- v) 
    return u_prima, v_prima


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

        u += u_prima * paso_h       # Recordar que u_prima = f
        v += v_prima * paso_h       # Recordar que v_prima = g

    return tiempo, valores_u, valores_v


## Trato de resolver el problema:

y0 = 0  # Condición inicial para y
v0 = 0  # Condición inicial para v


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



# ----------------------------- EULER IMPLICITO --------------------------------

def paso_siguiente(u,v,paso):
    """
    Calculo el paso siguiente de Un y Vn, es decir, Un+1 y Vn+1
    """
    v_sig = (v+(paso*k*(c-u)+lambda_*c_prima)*(1/m))/(1+(lambda_+(paso*paso)*k)*(1/m))
    u_sig = u + paso*v_sig

    return u_sig , v_sig


def euler_implicito(funcion, u0, v0, paso_h, intervalo_t):
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

        u, v = funcion(u, v, paso_h)


    return tiempo, valores_u, valores_v




# ----------------------------- RUNGE KUTTA --------------------------------

def calculo_q(u,v,paso):
    """
    Calculamos las q's
    """
    q1u = paso*v
    q1v = paso*((1/m)*(k*(c-u)+lambda_*(c_prima-v)))

    q2u = paso*(v+q1v)
    q2v = paso*((1/m)*(k*(c-(u+q1u))+lambda_*(c_prima-(v+q1v))))

    return q1u,q2u,q1v,q2v

def runge_kutta(funcion, u0, v0, paso_h, intervalo_t):
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

        q1u,q2u,q1v,q2v = calculo_q(u,v,paso_h)

        u += 0.5*(q1u+q2u)
        v += 0.5*(q1v+q2v)

    return tiempo, valores_u, valores_v



# ----------------------------- ORDEN DE CONVERGENCIA --------------------------------

def modulo(x1):
    """
    Operacion modulo de un numero
    """
    return raiz((x1**2))

def error_iteracion(x_actual , x_siguiente):
    """
    Calcula el error de iteracion
    """
    return modulo(x_actual-x_siguiente)/modulo(x_siguiente)

def orden_convergencia(iteracion):
    """
    Devuelve una lista con la convergencia calculada para cada paso.
    """
    orden_cv = []
    for i in range(1,iteracion-1):
        cv_iter_i = np.log(error_iteracion((iteracion+1)/iteracion))/ np.log((iteracion)/(iteracion-1))
        orden_cv.append(cv_iter_i)
    
    return orden_cv









