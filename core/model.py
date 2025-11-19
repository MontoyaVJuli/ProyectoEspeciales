import numpy as np
from scipy.sparse import spdiags, linalg
from core.config import ModelParameters, TreatmentParameters

class TumorModel:

    def __init__(self, params: ModelParameters):
        self.params = params
        self.N = params.grid_size
        self.grid = self._initialize_grid()
        self.tumor_mass_history = []

        self.r = self.params.D * self.params.dt / (2 * self.params.dx ** 2)
        self._setup_adi_matrices()

    def _initialize_grid(self):
        grid = np.zeros((self.params.grid_size, self.params.grid_size))
        center = self.params.grid_size // 2
        radius = self.params.initial_radius
        concentration = self.params.initial_concentration

        for i in range(self.params.grid_size):
            for j in range(self.params.grid_size):
                if (i - center) ** 2 + (j - center) ** 2 < radius ** 2:
                    grid[i, j] = concentration
        return grid

    def _apply_treatment(self, t: float, treatment: TreatmentParameters) -> float:
        if t < treatment.start_day:
            return 0.0

        cycle_length = treatment.duration_days + treatment.days_between_cycles

        time_in_cycle = (t - treatment.start_day) % cycle_length
        current_cycle_num = (t - treatment.start_day) // cycle_length

        if current_cycle_num < treatment.cycles and time_in_cycle < treatment.duration_days:
            return treatment.dose

        return 0.0

    #  Métodos para el Solver ADI (Crank-Nicolson)
    def _setup_adi_matrices(self):
        #Crea y factoriza las matrices tridiagonales dispersas A y B.
        N = self.N
        r = self.r

        diagonals_A = [np.full(N, -r), np.full(N, 1 + 2 * r), np.full(N, -r)]
        A = spdiags(diagonals_A, [-1, 0, 1], N, N, format='csc')

        A[0, 1] = -2 * r
        A[N - 1, N - 2] = -2 * r

        self.solver_A = linalg.factorized(A)
        self.solver_B = linalg.factorized(A)

    def _rk4_step(self, C_in, t, dt, treatment):

        def f(C, t_current):
            treatment_dose = self._apply_treatment(t_current, treatment)
            growth = self.params.k * C * (1 - C / self.params.K)
            ablation = treatment_dose * C
            return growth - ablation

        # 4 etapas rk
        k1 = dt * f(C_in, t)
        k2 = dt * f(C_in + 0.5 * k1, t + 0.5 * dt)
        k3 = dt * f(C_in + 0.5 * k2, t + 0.5 * dt)
        k4 = dt * f(C_in + k3, t + dt)

        return C_in + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

    def step(self, t: float, treatment: TreatmentParameters):
        #Avanza la simulación un paso de tiempo (dt)

        N = self.N
        r = self.r

        C_rk = self._rk4_step(self.grid, t, self.params.dt, treatment)

        # 2. PASO DE DIFUSIÓN (ADI Implícito sobre C_rk)
        C_star = np.zeros_like(C_rk)
        C_new = np.zeros_like(C_rk)

        D_y_explicit = r * (np.roll(C_rk, 1, axis=0) - 2 * C_rk + np.roll(C_rk, -1, axis=0))
        D_y_explicit[0, :] = 2 * r * (C_rk[1, :] - C_rk[0, :])  # Neumann top
        D_y_explicit[N - 1, :] = 2 * r * (C_rk[N - 2, :] - C_rk[N - 1, :])  # Neumann bottom
        b_x = C_rk + D_y_explicit

        for j in range(N):
            C_star[:, j] = self.solver_A(b_x[:, j])

        D_x_explicit = r * (np.roll(C_star, 1, axis=1) - 2 * C_star + np.roll(C_star, -1, axis=1))
        D_x_explicit[:, 0] = 2 * r * (C_star[:, 1] - C_star[:, 0])  # Neumann left
        D_x_explicit[:, N - 1] = 2 * r * (C_star[:, N - 2] - C_star[:, N - 1])  # Neumann right
        b_y = C_star + D_x_explicit

        for i in range(N):
            C_new[i, :] = self.solver_B(b_y[i, :])

        self.grid = np.clip(C_new, 0.0, self.params.K)

        # Guardar historial de frames
        self.tumor_mass_history.append(np.sum(self.grid))

    def run_simulation(self, total_days: int, treatment: TreatmentParameters):
        #Ejecuta la simulación completa
        self.tumor_mass_history = []
        steps_per_day = int(1 / self.params.dt)
        total_steps = total_days * steps_per_day

        frames = []

        for step in range(total_steps):
            current_day = (step * self.params.dt)
            self.step(current_day, treatment)

            # Guardar un "frame" cada día para la visualización
            if step % (steps_per_day * 1) == 0:
                frames.append(self.grid.copy())

        return frames, self.tumor_mass_history