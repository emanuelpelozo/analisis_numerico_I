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

y0 = 0  # Condición inicial para y
v0 = 0  # Condicion inicial para v