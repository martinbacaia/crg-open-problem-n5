# Fase 2 — Paso 30: derivación SIMBÓLICA EXACTA (sympy) de la cascada de cancelación Regge en
# grado 9 y grado 11 — piso exacto probado en grado 9 (s³), y un hallazgo mayor en grado 11
# (una familia de dimensión ≥8-9 ya satisface el bound s²), verificado en 2 configuraciones
# independientes de parámetros

## Motivación

Las partes 27 y 29 dejaron una cascada de cancelación numérica (`s⁵→s⁴→~s^3.5` en grado 11,
`s⁴→s³→?` en grado 9) que no lograba cerrar con precisión numérica suficiente (Richardson/mpmath se
quedaba sin resolución en el segundo o tercer escalón). Esta parte implementa lo que la parte 29
recomendaba explícitamente: una derivación **simbólica exacta** de los coeficientes de cancelación,
en vez de seguir extrapolando numéricamente.

## Método: por qué la cinemática de Regge de n=5 SÍ se puede tratar exactamente en sympy

Un primer intento directo (sustituir `E=1/x` en `regge_limit_n5.py` y llamar `sympy.series()` sobre
las expresiones tal cual) resultó completamente intratable: expresiones con `sqrt()` anidados dentro
de fracciones con más `sqrt()` (de la fórmula de boost de Lorentz) hacen que `series()` tarde minutos
u horas por término, y hay que evaluar cientos de términos (120 permutaciones × varios invariantes).

**Observación clave que lo hace tratable**: en toda la construcción de `build_kinematics_n5`
(`p1..p5`), la única fuente de irracionalidad algebraica no trivial es `sin_theta` (llamado `S` acá)
— y `S` aparece **siempre linealmente** (nunca elevado al cuadrado sin que ese cuadrado se pueda
sustituir de inmediato por `R(x) := 1-cos_theta(x)²`, una función racional conocida de `x`) en cada
componente de cada `p_i`. Esto se verifica explícitamente rastreando la construcción: `p3` tiene un
único componente con `S` (lineal); `P45 = -(p1+p2+p3)` hereda un único componente lineal en `S`; el
boost de Lorentz de `p4_rest,p5_rest` (que también tienen un único componente lineal en `S`, vía
`n45`) sólo produce productos de a lo sumo dos factores lineales en `S`, que se reducen de inmediato
a `S²=R(x)` (racional). Consecuencia: **toda la cinemática vive exactamente en la extensión
cuadrática `Q(x)(S)`, con `S²=R(x)`** — cerrada bajo `+,-,×,÷` (la división se resuelve
racionalizando por el conjugado, `(A+BS)/(C+DS) = (A+BS)(C-DS)/(C²-D²R)`, un truco algebraico
estándar). Se implementó esta aritmética exacta en `scripts/qext.py` (clase `QExt`, con `cancel()`
tras cada operación para evitar el crecimiento descontrolado de fracciones anidadas).

Con esto, `p_i`, `ε_i` (y por lo tanto `FdotF`, `Cabc`, `Msandwich`, y cualquier invariante) se
calculan como pares exactos `(A(x), B(x))` de funciones racionales de `x` — **sin ningún `sqrt()`
literal en ningún paso intermedio**. Sólo al final se sustituye `S` por su serie de potencias real
(`sympy.series(sqrt(R(x)), x, 0, N)`, rápido porque `R(x)` es una función racional simple, sin nada
anidado) para obtener la serie de Laurent física de cada invariante.

**Elección de parámetros**: se usan `t0`, `m0sq` como `sp.Rational` (en vez de floats) elegidos
perfectos cuadrados/simples (p.ej. `m0sq=4`), y el ángulo `phi0` (orientación fija del decaimiento
del sistema de retroceso `{4,5}`, ver parte 22) como un punto racional del círculo unitario
(`cos_phi0=3/5, sin_phi0=4/5`) en vez de un ángulo genérico — así **todo el resultado final queda en
`Q` (números racionales), sin ningún `sqrt` de constante suelto**, facilitando comparar resultados
entre configuraciones distintas.

**Optimización adicional decisiva**: en vez de recalcular `eta(a,b)` desde cero para cada una de las
120 permutaciones × cada invariante, se cachean **una sola vez** todos los `eta(p_i,p_j)`,
`eta(p_i,ε_j)`, `eta(ε_i,ε_j)` (25+25+25 valores), y a partir de ahí todos los `F^i:F^j` (10 pares),
`C^k_{ab}` (60 tríos), `M^{ij}_{ab}` (120 cuádruplas) — cada uno también una sola vez. Esto bajó el
tiempo de un cálculo de grado 11 completo (15 invariantes × 120 permutaciones) de "intratable" a
~6-9 minutos.

## Verificación cruzada obligatoria (antes de confiar en cualquier resultado)

1. **Cinemática**: se evaluó la construcción simbólica (`QExt` con `S` sustituido por su valor
   `sqrt(R(x))` exacto, no la serie) en `x=1/1000` y se comparó componente a componente contra
   `regge_limit_n5.py::build_kinematics_n5` (numérico, ya validado en la parte 23). Coincide a
   precisión de máquina (`~1e-11` a `~1e-14` de error relativo) para `p1..p5`, usando el mismo
   `phi0=0.3` que el script numérico (al principio se probó con `phi0` racional distinto y dio
   discrepancia — diagnosticado correctamente como una diferencia real de parámetro, no un bug, al
   volver a evaluar con el mismo `phi0=0.3` numérico).
2. **Polarizaciones**: mismo chequeo para `ε_1..ε_5` (construidas vía la misma proyección de gauge
   `growing_partner`), con semillas fijas — coincide a `~1e-13`-`~1e-15`.
3. Con esto verificado, se procedió a construir los invariantes de grado 9 y 11 con confianza en la
   cinemática base.

## Resultado grado 9 (B1, B2, B3 — dimensión 3, parte 27)

Se calcularon `B1(x), B2(x), B3(x)` como series de Laurent **exactas** (coeficientes racionales
exactos, no truncamiento numérico) sumando sobre las 120 permutaciones. Ejemplo (config. de
referencia `t0=-1, m0sq=4, cos_phi0=3/5, sin_phi0=4/5`):

```
B1: x^-8 = 82944/125    x^-7 = -1533952/125    x^-6 = 31146624/625   ...
B2: x^-8 = 41472/125    x^-7 = -766976/125     x^-6 = -2956608/625   ...
B3: x^-8 = -205568/125  x^-7 = 397888/125      x^-6 = -3237856/625   ...
```

**A1/A2 = 82944/41472 = 2 EXACTO** (confirma de forma exacta, no numérica, el hallazgo de la parte
27). **Hallazgo nuevo**: el cociente de los coeficientes en el orden SIGUIENTE (`x^-7`) es también
EXACTAMENTE 2 (`-1533952/-766976 = 2`) — es decir, `B1 - 2·B2` cancela **automáticamente** tanto el
término líder (`x^-8`) como el subsiguiente (`x^-7`), sin que eso se haya impuesto a mano.

**Análisis riguroso vía rango/núcleo exacto** (en vez de "probar una combinación y ver"): se armó la
matriz de coeficientes `3×5` (filas = órdenes `x^-8..x^-4`, columnas = `B1,B2,B3`) y se calculó el
rango exacto (sympy, `Matrix.rank()`, aritmética racional) acumulando filas en orden:

| orden incluido hasta | rango acumulado | dim. núcleo (combinaciones que cancelan hasta ahí) |
|---|---|---|
| `x^-8` | 1 | 2 |
| `x^-7` | 2 | 1 |
| `x^-6` | 3 | **0** |

El núcleo llega a dimensión 0 exactamente en `x^-6` — es decir, **ninguna combinación lineal no nula
de B1,B2,B3 puede anular también el término `x^-6`**, y la única dirección (hasta escala) que
sobrevive tras cancelar `x^-8` y `x^-7` es precisamente `B1-2B2`, cuyo primer término no nulo es
`x^-6` (=`E⁶`=`s³`). **Esto PRUEBA de forma exacta (no sugiere, no extrapola) que `s³` es el piso
real de la totalidad del espacio de grado 9** (que tiene dimensión exacta 3, parte 27 — B1,B2,B3 lo
agotan) — no sólo "el mejor resultado encontrado", sino el óptimo demostrado.

**Verificación de independencia de parámetros**: se repitió todo el cálculo con una configuración de
parámetros completamente distinta (`t0=-4, m0sq=9, cos_phi0=5/13, sin_phi0=12/13`, semillas de gauge
distintas). Resultado: **A1/A2=2 exacto de nuevo**, y la tabla de rango/núcleo da **exactamente el
mismo patrón** (núcleo 2→1→0 en `x^-8,x^-7,x^-6`) — el piso `s³` es robusto, no un artefacto de la
configuración numérica particular. (El coeficiente de la segunda dirección de cancelación, `c₃` en
`B2+c₃B3`, sí depende de los parámetros — `162/803` en la config. de referencia vs. `887/823` en la
segunda — pero eso es esperable: es la solución de una ecuación lineal con coeficientes que dependen
de la cinemática, no una identidad, y de todos modos no mejora sobre `B1-2B2`, según el propio
análisis de rango.)

## Resultado grado 11 (S1..S15 — dimensión ≥15, parte 28) — hallazgo mayor

Se repitió el mismo programa para los 15 invariantes base de la parte 28. Confirmado de forma EXACTA
(no sólo "6+ cifras" como en la parte 29):

- **A(S8) = -1/2 · A(S1)** exacto.
- **A(S14) = A(S15)** exacto.

(ambas relaciones, además, se re-confirmaron IDÉNTICAS en la segunda configuración de parámetros).

**Análisis de rango/núcleo exacto** (matriz `15×16`, filas = órdenes `x^-10` hasta `x^5`, columnas =
`S1..S15`), acumulando fila por fila:

| orden incluido hasta | rango acumulado | dim. núcleo |
|---|---|---|
| `x^-10` | 1 | 14 |
| `x^-9` | 2 | 13 |
| `x^-8` | 3 | 12 |
| `x^-7` | 4 | 11 |
| `x^-6` | 5 | 10 |
| `x^-5` | 6 | **9** |
| `x^-4` | 7 | 8 |
| `x^-3` | 8 | 7 |
| `x^-2` | 9 | 6 |
| `x^-1` | 10 | 5 |
| `x^0`  | 11 | 4 |
| `x^1`  | 12 | 3 |
| `x^2`  | 13 | 2 |
| `x^3`  | 14 | 1 |
| `x^4`  | 15 | **0** |

El rango sube exactamente en 1 en cada uno de los 15 pasos (sin degeneraciones extra más allá de las
2 relaciones exactas ya notadas, que son relaciones puntuales entre coeficientes individuales, no
relaciones de rango de fila). Lectura física: **el subespacio de combinaciones cuyos primeros 6
coeficientes (`x^-10` a `x^-5`) se anulan tiene dimensión 9** — es decir, existe una familia de
**dimensión 9** dentro del espacio completo de grado 11 (dimensión ≥15) cuyo crecimiento en este
límite de Regge es `≤ E⁴ = s²`, **satisfaciendo el bound conjeturado de Classical Regge Growth**, en
vivo contraste con la conclusión anterior (parte 24, ya marcada como provisional en la parte 28) de
que grado 11 topaba en `s⁴` — esa conclusión estaba basada en una familia de sólo 3 direcciones
(previa a incluir sistemáticamente el sandwich), muy por debajo de la dimensión real.

**Esta tabla se reprodujo IDÉNTICA (mismo rango en cada orden, mismo `x^4` como punto donde el
núcleo llega a 0) en la segunda configuración de parámetros** (`t0=-4, m0sq=9, cos_phi0=5/13,
sin_phi0=12/13`, semillas distintas) — evidencia fuerte de que no es un accidente numérico de una
configuración particular.

## Estado honesto — qué prueba esto y qué NO

**Lo que está probado con rigor (exacto, dos configuraciones independientes)**:
- Grado 9: el piso de crecimiento de TODO el espacio (dimensión 3) es exactamente `s³`, no hay
  ninguna combinación que baje más — demostrado, no sólo "no encontrado nada mejor".
- Grado 11: existe una familia de dimensión (al menos) 9 dentro del espacio de dimensión ≥15 cuyo
  crecimiento en ESTA trayectoria concreta de Regge es `≤ s²` — un resultado mucho más favorable que
  todo lo encontrado antes en el proyecto.

**Lo que NO está probado todavía (limitaciones reales, no retórica)**:
1. **Una sola familia de trayectorias de Regge**: el límite usado (parte 22/23) es UNA
   generalización concreta y ya documentada como "elección de modelado propia, no está en la
   literatura" del límite 2→2 de n=4 a n=5 (cluster duro `{1,2,3}` + retroceso `{4,5}`, con `t₁₃` y
   `M₄₅²` fijos). El bound CRG debe cumplirse para TODA generalización razonable de "s→∞ a
   invariantes de tipo t fijos" — round-tripar sólo 2 configuraciones numéricas de esta misma familia
   (aunque con parámetros `t0,m0sq,phi0` distintos) NO cubre otras posibles familias de límites de
   Regge para n=5 (p.ej. otras particiones en clusters). Esto es evidencia fuerte, no una prueba
   completa de que la familia de 9 dimensiones satisface CRG en general.
2. **Dimensión real de grado 11 podría ser >15** (parte 28 ya lo advertía: familias C y D usadas no
   exhaustivas). Si hay más generadores independientes, la familia de "grado ≤s²" podría ser aún más
   grande — o, en el peor caso, algún generador adicional podría acoplarse y estropear alguna
   combinación de la familia ya encontrada (aunque esto último es poco probable, ya que agregar más
   columnas a una matriz de rango completo no puede REDUCIR el rango de las columnas ya incluidas).
3. **No se verificó** si los elementos de esta familia de dimensión 9 son combinaciones
   "sensibles"/degeneradas (coeficientes gigantes que se cancelan de forma delicada) — sólo se probó
   que existen exactamente, en aritmética racional exacta, sin aproximación.
4. **No se exploró** si dentro de esta familia de dimensión 9 hay margen para bajar aún más allá de
   `s²` (el análisis se detuvo en el punto donde el núcleo llega a dimensión 0, que ocurre en
   `x^4` — de hecho el núcleo en `x^-4` ya es 8, es decir después de `x^-4` inclusive quedan 8
   direcciones con crecimiento `<s²` estrictamente, varias de las cuales decaen bien por debajo del
   bound en esta trayectoria particular — un resultado que en sí mismo es sorprendente y merece
   escrutinio antes de anunciarlo como logro central, dado el precedente de esta sesión de encontrar
   errores propios en resultados "demasiado buenos").

## Próximo paso concreto (por prioridad)

1. **Identificar explícitamente una base concreta (9 vectores) de la familia de grado 11 que satisface
   `≤s²`** en esta trayectoria (hasta ahora sólo se probó la DIMENSIÓN vía rango, no se extrajeron los
   vectores) y evaluarlos en una FAMILIA DISTINTA de límites de Regge (otro reparto de clusters, p.ej.
   `{1,2}` vs `{3,4,5}` con distinta partición interna) para ver si el bound se sostiene ahí también
   — el punto (1) de la sección anterior es la prioridad real antes de anunciar nada.
2. Ampliar la búsqueda de generadores de grado 11 más allá de las familias A+B+C+D (parte 28) para
   acotar mejor la dimensión real del espacio completo.
3. Si el punto 1 se sostiene en al menos una segunda familia de límites de Regge independiente,
   sería la primera evidencia seria de un candidato de grado 11 compatible con CRG para n=5 con spin
   — el objetivo original del proyecto.

## Archivos

- `scripts/qext.py` — aritmética exacta en la extensión cuadrática `Q(x)(S)`, `S²=R(x)`.
- `scripts/regge_n5_simbolico_kinematics.py` — cinemática de Regge n=5 exacta vía `QExt` (con
  comentario inline sobre cómo recalibrar parámetros para la verificación cruzada).
- `scripts/regge_n5_simbolico_eps.py` — polarizaciones `ε_i` exactas (gauge-fijadas) vía `QExt`.
- `scripts/regge_n5_simbolico_grado9.py` — cachea `eta`, `F:F`, `C`, `M` una sola vez; calcula
  `B1,B2,B3` como series de Laurent exactas (sumando las 120 permutaciones).
- `scripts/regge_n5_simbolico_grado11.py` — igual, para `S1..S15`.
- El análisis de rango/núcleo (tablas de este documento) se hizo con scripts ad hoc de pocas líneas
  (sympy `Matrix.rank()` sobre los coeficientes extraídos de los `.pkl` que generan los scripts de
  arriba) — no se dejaron como archivo separado; reproducibles en 2 minutos a partir de los `.pkl`
  (`B1,B2,B3` grado 9, `S1..S15` grado 11) que producen los scripts de arriba al correrlos.
