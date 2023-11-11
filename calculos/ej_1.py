import numpy as np
import matplotlib.pyplot as plt

# Parámetros dados
m = 1000  # masa en kg
k = 25000  # constante elástica en N/m
lambda_ = 0  # constante de amortiguación en Ns/m
c = 0.1  # elevación del terreno en metros

# Condiciones iniciales
y0 = 0.1  # posición inicial
v0 = 0  # velocidad inicial

# Tiempo
t_final = 5  # segundos