import streamlit as st
from core.model import TumorModel, ModelParameters, TreatmentParameters
from utils.visualization import plot_tumor_grid, plot_growth_curve

st.set_page_config(layout="wide", page_title="Simulador Terapia Oncológica")

st.title(" Simulador De Evolución Tumoral Según Terapia Oncológica")
st.markdown(
    "Simulador para evaluar el efecto de terapias oncológicas  _(Quimio/Radioterapia)_  en el crecimiento de las células tumorales basado  en  _**Ecuaciones diferenciales Parciales (EDPs) de Reacción-Difusión**_  con un  **Solver Semi-Implícito**  _(Crank-Nicolson para Difusión y Runge-Kutta 4 para Reacción)_ de alta estabilidad y precisión.")

#  Inicializar Session State para guardar resultados
if 'simulation_run' not in st.session_state:
    st.session_state.simulation_run = False
if 'frames' not in st.session_state:
    st.session_state.frames = []
if 'history' not in st.session_state:
    st.session_state.history = []
if 'sim_days_total' not in st.session_state:
    st.session_state.sim_days_total = 0
if 'model_params' not in st.session_state:
    st.session_state.model_params = ModelParameters()

#simulación con caché
@st.cache_data(show_spinner="Ejecutando simulación. Esto puede tardar unos segundos...")
def run_cached_simulation(model_params: ModelParameters, treatment_params: TreatmentParameters, sim_days: int):
    model = TumorModel(model_params)
    frames, history = model.run_simulation(sim_days, treatment_params)
    return frames, history

# Contenedor Principal
tab_config_tumor, tab_config_treat, tab_results = st.tabs([
    "Propiedades Biológicas del Tumor",
    "Configuración del Tratamiento _(Protocolo)_",
    "Simulación de Resultados"
])

# Configuración del Tumor
with tab_config_tumor:
    st.header("Configuración de Parámetros del Tumor   _(componentes biológicos)_")
    st.info("Ajuste los parámetros del tumor para modelar diferentes tipos de cáncer.")

    col_tumor_1, col_tumor_2 = st.columns(2)

    with col_tumor_1:
        sim_days = st.slider("Duración de la Simulación (Días)", 50, 500, 200,
                             help="Tiempo total del experimento simulado.")

        initial_radius = st.slider("Radio Inicial del Tumor ", 2, 20, 5,
                                   help="Tamaño inicial del tumor en el día 0.")

        initial_concentration = st.slider("Concentración Inicial", 0.1, 1.0, 0.8, format="%.2f",
                                          help="Densidad celular inicial en el núcleo del tumor (1.0 = Máxima).")

    with col_tumor_2:
        D_coeff = st.slider("Velocidad de Esparcimiento", 0.01, 0.2, 0.05, format="%.2f",
                            help="Capacidad de las células tumorales para migrar e invadir el tejido sano adyacente.")

        k_growth = st.slider("Velocidad de Crecimiento", 0.005, 0.05, 0.02, format="%.3f",
                             help="Velocidad con que las células se reproducen. Mayor Velocidad = tumor más agresivo.")
        grid_size = 100

# Configuración del Tratamiento
with tab_config_treat:
    st.header("Esquema de Quimio/Radioterapia")
    st.info("Defina el protocolo de dosificación y su efectividad.")

    col_treat_1, col_treat_2 = st.columns(2)

    with col_treat_1:
        treat_start = st.slider("Día de Inicio de Terapia", 0, 100, 50, help="Momento en que comienza la intervención.")
        treat_dose = st.slider("Efectividad / Dosis Máxima", 0.0, 0.1, 0.05, format="%.2f", help="Potencia del tratamiento para inducir muerte celular.")
        treat_cycles = st.slider("Nº de Ciclos", 0, 10, 3, help="Número de veces que se repite el tratamiento.")


    with col_treat_2:
        treat_duration = st.slider("Duración del Ciclo _(días)_", 1, 10, 5, help="Número de días que dura un ciclo del tratamiento.")
        treat_pause = st.slider("Días de Recuperación ", 5, 20, 10, help="Número de días de pausa entre los ciclos del tratamiento.")

#  TAB 3:Resultados
with tab_results:
    st.header("Simulación de Resultados")

    # 1. Botón de Ejecución (Guarda en el Session State)
    if st.button(" >️ Iniciar Simulación ", type="primary"):
        current_model_params = ModelParameters(
            grid_size=grid_size,
            D=D_coeff,
            k=k_growth,
            initial_radius=initial_radius,
            initial_concentration=initial_concentration
        )

        current_treatment_params = TreatmentParameters(
            start_day=treat_start,
            dose=treat_dose,
            cycles=treat_cycles,
            duration_days=treat_duration,
            days_between_cycles=treat_pause
        )

        #Ejecuta simulación cache.
        frames, history = run_cached_simulation(current_model_params, current_treatment_params, sim_days)

        #Guardar resultados en el estado de la sesión
        st.session_state.simulation_run = True
        st.session_state.frames = frames
        st.session_state.history = history
        st.session_state.sim_days_total = sim_days
        st.session_state.model_params = current_model_params

        st.success("_Simulación completada,_ navegue los resultados a continuación.")

    # 2. Visualización de Resultados
    if st.session_state.simulation_run:

        frames = st.session_state.frames
        history = st.session_state.history
        sim_days_total = st.session_state.sim_days_total
        model_params = st.session_state.model_params

        initial_mass = history[0]
        final_mass = history[-1]

        delta_mass = final_mass - initial_mass
        change_percentage = 100 * delta_mass / initial_mass if initial_mass > 0 else 0

        if delta_mass < 0:
            delta_text = f"- {abs(change_percentage):.2f}% Reducción"
        else:
            delta_text = f"+ {change_percentage:.2f}% Crecimiento"

        st.subheader("Resumen del Tratamiento")
        metric_col_1, metric_col_2, metric_col_3 = st.columns(3)

        metric_col_1.metric("_Masa Inicial (Día 0)_", f"{initial_mass:.2f} u.a.")

        metric_col_2.metric(
            label="Masa Final",
            value=f"{final_mass:.2f} u.a.",
            delta=delta_text,
            delta_color="inverse"
        )

        metric_col_3.metric("Días de Terapia Aplicada", f"{treat_cycles * treat_duration} días")


        st.subheader(" Vista Inicial y Final del Tumor")
        col_static_1, col_static_2 = st.columns(2)

        with col_static_1:
            st.markdown("###### Inicio de la Simulación Día 0")
            fig_initial = plot_tumor_grid(frames[0], day=0)
            st.pyplot(fig_initial)

        with col_static_2:
            st.markdown(f"###### Final de la Simulación Día {sim_days_total}")
            fig_final = plot_tumor_grid(frames[-1], day=sim_days_total)
            st.pyplot(fig_final)

        st.subheader("Evolución Temporal")

        viz_col_1, viz_col_2 = st.columns(2)

        dt = model_params.dt
        frame_interval_days = int(1 / dt) * 1
        max_days = sim_days_total

        with viz_col_1:
            st.markdown("###### Vista del Tumor _(Mapa de Calor)_")

            selected_day = st.slider("Navegar por Días de Simulación:", 0, max_days, max_days, step=frame_interval_days)

            selected_frame_index = selected_day

            if selected_frame_index >= len(frames):
                selected_frame_index = len(frames) - 1

            fig_anim = plot_tumor_grid(frames[selected_frame_index], day=selected_day)
            st.pyplot(fig_anim)


        with viz_col_2:
            st.markdown("###### Curva de Evolución  _(Masa total en el tiempo)_")
            fig_curve = plot_growth_curve(history, dt)

            treat_start = treat_start
            fig_curve.axes[0].axvline(x=treat_start, color='red', linestyle='--', alpha=0.7, label='Inicio Terapia')
            fig_curve.axes[0].legend()

            st.pyplot(fig_curve)



    elif not st.session_state.simulation_run:
        st.info("Presiona el botón para ejecutar la simulación y ver los resultados.")