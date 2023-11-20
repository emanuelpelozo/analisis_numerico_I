import numpy as np
import matplotlib.pyplot as plt
from ej_1 import pasos_h, m, k, y0,v0, INTERVALO_T

nuevo_lambda = 750
lambda_ = nuevo_lambda
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




# ----------------------------- K Y LAMBDA OPTIMOS --------------------------------


"""
Vamos a buscar un K y un lambda, tal que el sistema de amortiguacion no se comprima 
mas de 0.05m, es decir y-c >= 0.05m. 
Para buscar el optimo lo haremos por Fuerza Bruta, primero con un lambda constante
y con distintos valores de k hasta encontrar el optimo k, para luego ir variando 
lambda y con el k optimo, hasta en contrar el optimo de lambda. Luego haremos
al reves, primero buscaremos el optimo lambda y luego el optimo k.
"""

def calculo_q(u,v,paso,k,lambda_, c, c_prima):
    """
    Calculamos las q's
    """
    q1u = paso*v
    q1v = paso*((1/m)*(k*(c-u)+lambda_*(c_prima-v)))

    q2u = paso*(v+q1v)
    q2v = paso*((1/m)*(k*(c-(u+q1u))+lambda_*(c_prima-(v+q1v))))

    return q1u,q2u,q1v,q2v


def optimo(max_compresion):
    #Vemos cuales valores cumplen (su compresion maxima es mayor que -0.05) y nos quedamos con el mas grande
    validos = []
    for compresion in max_compresion:
        if compresion[1] >= -0.05:
            validos.append(compresion)
    
    optimo = max(validos, key = lambda item:item[1])
    return optimo


def optimoLambdaComp(paso_h, intervalo_t, u0, v0, k_opt):
    """
    Como lambda tiene valores mas pequenios que K realizamos lo mismo pero
    con valores de a 150, de 150 hasta 3000.
    """
    valores_lambda = []
    for i in range(150,3150,150): valores_lambda.append(i)
    max_compresion = []

    #No queria hacerle modificaciones al RK original asi que copie y pegue el codigo para hacerle sus modificaciones aca.
    tiempo = np.arange(0, intervalo_t, paso_h)
    
    for lam in valores_lambda:
        minimo = float("inf")

        # seteamos condiciones iniciales
        u = u0
        v = v0

        # Recorremos los puntos en el tiempo.
        for t in tiempo:
            c , c_prima = valor_de_c(t)
            q1u,q2u,q1v,q2v = calculo_q(u,v,paso_h,k_opt,lam,c,c_prima)

            u += 0.5*(q1u+q2u)
            v += 0.5*(q1v+q2v)

            #Buscamos el valor de compresion mas grande (lo mas comprimido posible, es decir, cuando u esta mas cercano a c)
            if (u-c)<minimo:
                minimo = (u-c)

        max_compresion.append((lam,minimo))

    lambda_optimo = optimo(max_compresion)
    return lambda_optimo


def optimoLambdaAcel(paso_h, u0, v0, k_opt):
    """
    Como lambda tiene valores mas pequenios que K realizamos lo mismo pero
    con valores de a 150, de 150 hasta 3000.
    """
    valores_lambda = []
    for i in range(150,3150,150): valores_lambda.append(i)
    max_aceleracion = []

    #No queria hacerle modificaciones al RK original asi que copie y pegue el codigo para hacerle sus modificaciones aca.
    tiempo = np.arange(1.1, 1.4, paso_h)
    
    for lam in valores_lambda:
        max_ac = float("-inf")

        # seteamos condiciones iniciales
        u = u0
        v = v0

        # Recorremos los puntos en el tiempo.
        for t in tiempo:
            c , c_prima = valor_de_c(t)
            q1u,q2u,q1v,q2v = calculo_q(u,v,paso_h,k_opt,lam,c,c_prima)

            u += 0.5*(q1u+q2u)
            v += 0.5*(q1v+q2v)

            #Buscamos el valor de aceleracion mas chico al pasar por la loma de burro, es decir, en dicho intervalo de tiempo.
            aceleracion = q1v/paso_h
            if aceleracion>max_ac:
                max_ac = aceleracion

        max_aceleracion.append((lam,max_ac))

    lambda_optimo = min(max_aceleracion, key = lambda item:item[1])
    return lambda_optimo



def optimoK(paso_h, intervalo_t, u0, v0):
    """
    Podemos ver que K tiene valores grandes, por lo que tomaremos valores de
    a mil, desde mil hasta 30 mil. Por cada uno, veremos cual es la maxima 
    compresion y luego buscaremos la minima compresion del total.
    """
    lambda_ = 750
    valores_k = []
    for i in range(1000,31000,1000): valores_k.append(i)

    max_compresion = []
    max_aceleracion = []

    #No queria hacerle modificaciones al RK original asi que copie y pegue el codigo para hacerle sus modificaciones aca.
    tiempo = np.arange(0, intervalo_t, paso_h)
    
    for k in valores_k:
        min_comp = float("inf")
        max_ac =   float("-inf")

        # seteamos condiciones iniciales
        u = u0
        v = v0

        # Recorremos los puntos en el tiempo.
        for t in tiempo:
            c , c_prima = valor_de_c(t)
            q1u,q2u,q1v,q2v = calculo_q(u,v,paso_h,k,lambda_,c,c_prima)

            u += 0.5*(q1u+q2u)
            v += 0.5*(q1v+q2v)

            #Buscamos el valor de compresion mas grande (lo mas comprimido posible, es decir, cuando u esta mas cercano a c)
            if (u-c)<min_comp:
                min_comp = (u-c)

            #Buscamos el valor de aceleracion mas chico al pasar por la loma de burro, es decir, en dicho intervalo de tiempo.
            if (t>=1.1) and (t<=1.4):
                aceleracion = q1v/paso_h
                if aceleracion>max_ac:
                    max_ac = aceleracion

        max_compresion.append((k,min_comp))
        max_aceleracion.append((k,max_ac))
    
    k_opt_comp = optimo(max_compresion)
    k_opt_acel = min(max_aceleracion, key = lambda item:item[1]) #Buscamos la aceleracion mas chica para todas las aceleraciones que hay.

    lambda_opt_comp = optimoLambdaComp(paso_h, intervalo_t, u0, v0, k_opt_comp)
    lambda_opt_acel = optimoLambdaAcel(paso_h, u0, v0, k_opt_acel)

    return k_opt_comp, k_opt_acel, lambda_opt_comp, lambda_opt_acel



def optimosFB(paso_h, intervalo_t, u0, v0):
    """
    Calculamos el lambda y K optimos a partir de Fuerza Bruta.
    """

    valores_k = []
    for i in range(1000,91000,1000): valores_k.append(i)

    valores_lambda = []
    for i in range(150,5250,150): valores_lambda.append(i)

    min_compresion = []     #Lista con tuplas de la forma (k,lambda,compresion)
    max_aceleracion = []    #Lista con tuplas de la forma (k,lambda,aceleracion)

    tiempo = np.arange(0, intervalo_t, paso_h)

    for k in valores_k:

        for lam in valores_lambda:

            min_comp = float("inf")
            max_ac =   float("-inf")

            u = u0
            v = v0

            for t in tiempo:
                c , c_prima = valor_de_c(t)
                q1u,q2u,q1v,q2v = calculo_q(u,v,paso_h,k,lam,c,c_prima)
                
                u += 0.5*(q1u+q2u)
                v += 0.5*(q1v+q2v)

                #Buscamos el momento en el que el resorte este mas comprimido.
                if (u-c)<min_comp:
                    min_comp = (u-c)

                #Buscamos la aceleracion maxima al pasar por la loma de burro.
                if (t>=1.1) and (t<=1.4):
                    aceleracion = q1v/paso_h
                    if aceleracion>max_ac:
                        max_ac = aceleracion

            if (min_comp)>=-0.05:
                min_compresion.append((k,lam,min_comp))

            max_aceleracion.append((k,lam,max_ac))
    
    opt_k_comp, opt_lam_comp, _ = max(min_compresion, key=lambda item:item[2])
    opt_k_ac, opt_lam_ac, _ = min(max_aceleracion, key=lambda item:item[2])

    return opt_k_comp, opt_k_ac, opt_lam_comp, opt_lam_ac