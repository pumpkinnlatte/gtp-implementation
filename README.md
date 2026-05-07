# packing-functions
Implementacion de las funciones de empaquetamiento estudiadas
en Torres-Jimenez et al. (2015).
 
Lista de algoritmos implementados
 
| # | Algoritmo                  |
|---|----------------------------|
| 1 | DirectGTP                  | 
| 2 | Binomial                   |
| 3 | InverseGTP                 | 
| 4 | NextBinomial               | 
| 6 | DirectDP                   |
| 7 | DirectBP                   | 
| 8 | InverseBP                  | 
 
## Estructura
 
    include/    headers (arith.h, dp.h, gtp.h, bp.h, packing_types.h)
    src/        implementacion
    tests/      suite de pruebas
    examples/   programas demostrativos
    build/      artefactos 
## Uso
 
    make                                    # compila todo
    make test                               
    make demo                               # ejecuta examples/demo
    make cli ARGS="dp encode 2 3"           
    make debug                              
    make clean
 
### CLI interactiva
 
El programa `p-functions` permite ejecutar los polinomios por linea de comandos.
 
Ejemplos:
 
    ./build/p-functions dp encode 2 3          # vector DP -> alpha = 17
    ./build/p-functions dp decode 2 17         # alpha -> vector DP = (2,3)
    ./build/p-functions gtp encode 0 2 5       # vector GTP -> alpha
    ./build/p-functions gtp decode 3 8         # alpha -> vector GTP
    ./build/p-functions bp encode 1 2 1        # vector BP -> alpha = 21
    ./build/p-functions bp decode 3 21         # alpha -> vector BP = (1,2,1)
    ./build/p-functions dp table 6             # imprime tabla 6x6 de DP
    ./build/p-functions help                   # ayuda completa

## Limites del tipo entero

Se usa un entero uint64_t (0 a 2^64 - 1) por lo que por ahora hay limitaciones
con dimensionalidades altas
