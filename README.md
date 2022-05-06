# Solver ecuación de la energía
## Código del trabajo de la asignatura de Mecánica de Fluidos Computacional y Experimental

Escrito por:
- Elfidio Ángel Alonso Artal 
- Jorge Fernández Chinchón
- Carlos González Pérez
- Diego Jordá Espí
- Luis Navarro García

## Instalación
Se recomienda crear un entorno virtual, con Anaconda o cualquier otro método, para mantener las dependencias aisladas.
El código utiliza **Python 3.10.4**, pero debería funcionar en versiones cercanas.

Para instalar el código, clonar el repositorio: 
```bash
git clone https://github.com/LuisNavarroGarcia/Trabajo_CFD.git
```
O, si no se desea clonar el repositorio, descargar el *.zip* de Github: https://github.com/LuisNavarroGarcia/Trabajo_CFD/archive/refs/heads/master.zip y descomprimirlo en la carpeta que se desee.

Después es necesario acceder a la carpeta mediante:
```bash
cd Trabajo_CFD
```

Para instalar las dependencias:
```bash
pip install -r requirements.txt
```

Relizando esto, el software debería estar listo para su uso.

## Utilización
Para ejecutar el código, desde la carpeta *Trabajo_CFD*:
```bash
python launcher.py
```

Éste comando ejecutará el solver con la configuración por defecto que se encuentra en el archivo *default_config.txt*. Este archivo se puede modificar para cambiar la configuración

El archivo *default_config.txt* se puede duplicar para crear distintos archivos de configuración. Como los archivos customizados se guardarán con otros nombres, si se desea utilizar un archivo de configuración customizado, se puede hacer con el comando:
```bash
python launcher.py --config_file [nombre del archivo con la extensión]
```

Ejemplo:
```bash
python launcher.py --config_file custom_config.txt
```

Donde *custom_config.txt* es el archivo de configuración creado por el usuario.

A consinuación se muestra el archivo *default_config.txt*, con la configuración por defecto:

    [fluid properties]
    k = 0.024
    cv = 718.5
    rho = 1.225
    
    [mesh]
    num_cells = 4
    num_bc = 4
    
    [type boundary conditions]
    boundary_1 = neumann
    boundary_2 = neumann
    boundary_3 = neumann
    boundary_4 = neumann
    
    [value boundary conditions]
    boundary_1 = lambda x, y, t : 293
    boundary_2 = lambda x, y, t : 700
    boundary_3 = lambda x, y, t : 1100
    boundary_4 = lambda x, y, t : 450
    
    [initial conditions]
    u = lambda x, y, t: np.full((len(x), 2), np.array([0., 0.006]))
    t0 = 0
    
    [propagators]
    propagator = euler_implicit
    
    [integrators]
    diffusion_integrator = difusion_cds_weighted
    convection_integrator = conv_upwind_order1
    
    [time configuration]
    dt_calc = dt_constant
    dt0 = 0.01
    courant = 10
    
    [stopping criteria]
    stop_at_time_flag = 0
    stop_at_iteration_flag = 0
    stop_at_max_error_flag = 1
    stop_at_mean_error_flag = 1
    
    stop_at_time = 10
    stop_at_iteration = 1000
    stop_at_max_error = 0.007
    evaluate_max_error_at_n_last_iterations = 2
    stop_at_mean_error = 0.005
    evaluate_mean_error_at_n_last_iterations = 2
    
    OR_AND_condition_for_time_condition = 0
    OR_AND_condition_for_iteration_condition = 0
    OR_AND_condition_for_max_error_condition = 0
    OR_AND_condition_for_mean_error_condition = 0
    
    
    [plot configuration]
    active_convergence_plot = 1
    plot_max_error = 1
    plot_mean_error = 1
    plot_frecuency = 5
    
    contour_colormap = inferno
    contour_interpolation_method = cubic

En [fluid properties] se introducen:
- **k**: Conductividad térmica *[W/(m·K)]*
- **cv**: Calor específico *[J/(kg·K)]*
- **rho**: Densidad *[kg/m^3]*

En [mesh], se introducen:
- **num_cells**: Número de celdas de la malla
- **num_bc**: Número de condiciones de contorno

Los archivos de la malla y condiciones de contorno deben encontrarse dentro de la carpeta *mesh* y el formato de normbre de los archivos debe ser:
- malla: *cells_[# de celdas].dat*, y *nodes_[# de celdas].dat*.
- Condiciones de contorno: *bc_[# de condición]_[# de celdas]*

En [type of boundary condition] se introduce:
- **boundary_1**: neumann o dirichlet
- **boundary_2**: neumann o dirichlet
- **boundary_3**: neumann o dirichlet
- **boundary_4**: neumann o dirichlet
- **...**

En [value of boundary condition] se introduce:
- **boundary_1 = lambda x, y, t:** Función dependiente de *x*, *y*, *t*
- **boundary_2 = lambda x, y, t:** Función dependiente de *x*, *y*, *t*
- **boundary_3 = lambda x, y, t:** Función dependiente de *x*, *y*, *t*
- **boundary_4 = lambda x, y, t:** Función dependiente de *x*, *y*, *t*
- **...**

En [initial conditions] se introduce:
- **u = lambda x, y, t: np.full((len(x), 2), np.array([{1}, {2}]))**: Donde {1} magnitud de la velocidad en x, puede depender de *x*, *y*, *t*. Y {2} magnitud de la velocidad en y, puede depender de *x*, *y*, *t*.
- **t0** : Tiempo inicial

En [propagators] se introduce:
- **propagator**: Puede ser *euler_explicit*, *euler_implicit*, *euler_pred_corr*, *crank_nicolson* y *runge_kutta4*.

En [integrators] se introduce:
- **diffusion_integrator**: Puede ser *diffusion_cds* o *difusion_cds_weighted*
- **convection_integrator**: Puede ser *conv_upwind_order1* o *conv_cds*

En [time_configuration] se introduce:
- **dt_calc**: Puede ser *dt_constante* o *dt_adaptative*
- **dt0**: *dt* inicial
- **courant**: Constante de Courant

En [stopping criteria] se introduce:
- **stop_at_time_flag**: 1 para criterio de parada de tiempo activo, 0 inactivo.
- **stop_at_iteration_flag**: 1 para criterio de parada de iteración máxima activo, 0 inactivo.
- **stop_at_max_error_flag**: 1 para criterio de parada de error máximo, 0 inactivo.
- **stop_at_mean_error_flag**: 1 para criterio de parada de error medio activo, 0 inactivo.
- **stop_at_time**: Tiempo para parar *solver*
- **stop_at_iteration**: Iteración para parar *solver*
- **stop_at_max_error**: Error máximo para parar *solver*
- **evaluate_max_error_at_n_last_iterations**: Celdas donde se evalúa el error máximo (comenzando desde atrás).
- **stop_at_mean_error**: Error medio para parar *solver*
- **evaluate_mean_error_at_n_last_iterations**: Celdas donde se evalúa el error medio (comenzando desde atrás).
- 
Las siguientes *flags* indican si la condición que indica su nombre actúa un **OR** o un **AND**. Si es todo 0, el primer criterio que converja parará el *solver*, cuando alguna es 1, esa condición se tendrá que cumplir sí o sí.

- **OR_AND_condition_for_time_condition**: 0/1
- **OR_AND_condition_for_iteration_condition**: 0/1
- **OR_AND_condition_for_max_error_condition**: 0/1
- **OR_AND_condition_for_mean_error_condition**: 0/1

En [plot configuration] se introduce:
- **active_convergence_plot**: 1 para graficar la evolución del error para observar si se llega a la convergencia (más lento de ejecutar). 0 si no se desea graficar.
- **plot_max_error**: 1 si se plotea el error máximo, 0 si no.
- **plot_mean_error**: 1 si se plotea el error medio, 0 si no.
- **plot_frecuency**: Cada cuántas iteraciones se grafica la gráfica anteriormente mencionada.
- **contour_colormap**: Mapa de color para la gráfica de contornos. Puede ser *inferno*, *binary*, *viridis*, *plasma*, *magma* y *cividis*.
- **contour_interpolation_method**: Tipo de interpolación para gráfica de contorno. Puede ser: *cubic*, *nearest* y *linear*.