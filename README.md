# Simulador de Colisión de Fluidos

Este repositorio contiene una implementación experimental de un simulador de partículas para fluidos, orientado al estudio de colisiones y transferencias de masa en un volumen tridimensional. El código está escrito en C++20 y utiliza CMake como sistema de construcción.

## Características principales

- Lectura de archivos binarios `.fld` con las condiciones iniciales de la simulación (densidad de partículas, número de partículas y parámetros geométricos). 
- Construcción de una malla tridimensional de bloques (`Grid`) que agrupa partículas por posición para acelerar la búsqueda de vecinos y el cálculo de colisiones.
- Estructuras de datos para partículas (`Particle`) y vectores tridimensionales (`Vect3`) con operadores aritméticos básicos y utilidades para evaluar distancias.
- Parámetros físicos configurados en `sim/libraries.hpp`, como gravedad, viscosidad, tamaño de partícula o constantes de amortiguamiento empleadas en los cálculos.
- Rutinas para cargar y contrastar trazas binarias (`.trz`) de estados intermedios de la simulación, lo que permite validar los resultados frente a referencias externas.

## Estructura del proyecto

```
Fluid-Collision-Simulation/
├── CMakeLists.txt        # Configuración principal de CMake y dependencias (GoogleTest, GSL)
├── fluid/                # Ejecutable principal de la simulación
├── sim/                  # Biblioteca con la lógica de partículas, malla y utilidades
├── utest/                # Esqueleto de pruebas unitarias con GoogleTest
├── input.fld             # Ejemplo de archivo de entrada en formato binario
├── small.fld             # Otro ejemplo de entrada de menor tamaño
└── func_antiguas.txt     # Notas y funciones previas (histórico)
```

Cada subdirectorio posee su propio `CMakeLists.txt`, lo que permite compilar librerías y ejecutables de forma modular.

## Requisitos previos

- CMake ≥ 3.21
- Compilador compatible con C++20 (GCC 11+, Clang 13+, MSVC 19.3+)
- Git (para clonar el repositorio y permitir que CMake descargue GoogleTest y GSL)

No es necesario instalar GoogleTest ni GSL manualmente: CMake los descarga automáticamente durante la configuración.

## Compilación

1. Crear un directorio de construcción:
   ```bash
   cmake -S . -B build
   ```
2. Compilar los objetivos principales (`sim` y `fluid`):
   ```bash
   cmake --build build
   ```

Los binarios se ubican en `build/fluid/` (ejecutable `fluid`) y `build/sim/` (biblioteca estática/dinámica según la plataforma).

## Ejecución de la simulación

El ejecutable `fluid` espera tres argumentos:

```bash
./build/fluid/fluid <pasos> <archivo_entrada.fld> <prefijo_salida>
```

- `<pasos>`: entero que representa la cantidad de iteraciones o pasos de simulación solicitados.
- `<archivo_entrada.fld>`: ruta al archivo binario con las condiciones iniciales.
- `<prefijo_salida>`: prefijo para los archivos de salida (actualmente reservado para futuras extensiones).

Ejemplo utilizando el archivo de muestra incluido:

```bash
./build/fluid/fluid 1 input.fld resultados
```

Durante la ejecución, el programa informará por consola la lectura de parámetros (`ppm`, `np`) y la cantidad de partículas cargadas en memoria.

## Formato de los archivos `.fld`

El lector binario definido en `sim/progargs.cpp` espera la siguiente estructura:

1. Un número en coma flotante (`float`) con el valor de partículas por metro (`ppm`).
2. Un entero (`int`) con el número total de partículas (`np`).
3. Para cada partícula, nueve valores `float` consecutivos: posición `(x, y, z)`, pseudo-velocidad `hv` y velocidad `v`.

A partir de `ppm` y constantes globales (`const_r`, `global_density`), el programa deriva la masa y el tamaño de celda para crear la malla de bloques.

## Parámetros físicos y constantes

Todos los valores físicos que gobiernan la simulación (gravedad, viscosidad, rigidez de colisión, tamaño de partícula, etc.) se encuentran en `sim/libraries.hpp`. Modificar estas constantes permite ajustar el comportamiento del sistema sin alterar la lógica principal.

## Trazas binarias (`.trz`)

Las funciones disponibles en `sim/particles_motion.cpp` permiten cargar y validar archivos `.trz` que representan estados intermedios de la simulación (por ejemplo, densidad o aceleraciones transferidas). Estas trazas resultan útiles para comparar la salida del simulador con referencias generadas por otras herramientas.

## Pruebas unitarias

El directorio `utest/` contiene un esqueleto inicial de pruebas basado en GoogleTest. Para incluirlo en la compilación, añade el subdirectorio al flujo de CMake (`add_subdirectory(utest)`) y, posteriormente, ejecuta:

```bash
cmake --build build --target utest
ctest --test-dir build
```

Este conjunto de pruebas sirve como punto de partida para automatizar verificaciones sobre lectura de archivos, colisiones y dinámica de partículas.

## Próximos pasos sugeridos

- Completar y ampliar los casos de prueba en `utest/`.
- Implementar la escritura de resultados con el prefijo pasado como tercer argumento.
- Documentar formalmente el formato de archivos `.trz` y generar ejemplos de referencia.

## Licencia

El proyecto no especifica aún una licencia. Antes de distribuir o reutilizar el código, asegúrate de definir la licencia apropiada conforme a las necesidades del equipo.
