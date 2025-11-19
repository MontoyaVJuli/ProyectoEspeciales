from dataclasses import dataclass


@dataclass
class ModelParameters:
    grid_size: int = 100  # rejilla (N x N)
    D: float = 0.05  # coeficiente de difusion
    k: float = 0.02  # tasa de crecimiento (proliferacion)
    K: float = 1.0  # capacidad de carga maxima, lim biologico
    dt: float = 0.1  # paso de tiempo (Delta t)
    dx: float = 1.0  # paso espacial (Delta x, asumido 1.0)
    initial_radius: int = 5  # radio inicial del tumor
    initial_concentration: float = 0.8  # concentracion inicial en el nucleo del tumor


@dataclass
class TreatmentParameters:
    start_day: int = 20  # dia de inicio del tratamiento
    dose: float = 0.05  # efectividad/dosis de la quimio/radio
    duration_days: int = 5  # duracion de cada ciclo (días)
    cycles: int = 3  #  ciclos
    days_between_cycles: int = 10  # pausa entre ciclos (días)