# Fase 2 — Paso 33: backend rápido (lseries.py, sin sympy en el loop caliente) — bug real
# encontrado y corregido en el camino, verificación exacta de todo lo anterior, y el subespacio
# universal de dimensión 2 (parte 32) confirmado en una SEGUNDA configuración independiente

## Motivación

El usuario preguntó si convenía usar Singular u otro lenguaje más rápido que Python para acelerar
los cálculos de grado 11 (varios minutos por corrida). Evaluación: Singular no es la herramienta
correcta (su fuerza es álgebra de Gröbner/anillos de invariantes, no aritmética de series de Laurent
en 1 variable). Un benchmark directo mostró que el cuello de botella real es el overhead de árbol de
expresión de sympy para algo que es, en el fondo, convolución de listas de fracciones: multiplicar y
truncar una serie 2000 veces dio **sympy=106s vs. `fractions.Fraction` puro=1.6s (~65x)**.

## `lseries.py`: serie de Laurent truncada exacta sobre Q, con `fractions.Fraction`

Implementa `+,-,*,potencia entera` (sin división — ninguna fórmula de `eta`, `F:F`, `C`,
`Msandwich`, `S1..S15`/`B1..B3` la necesita). Verificado con `test_lseries.py` (100 casos aleatorios,
mezclando +,-,*,**2,**3, escalares, órdenes negativos y positivos) contra sympy como referencia.

**Bug real encontrado por el propio test** (v1 de `lseries.py`): truncar el orden superior después
de CADA multiplicación (para acotar tamaño) es INCORRECTO para series de Laurent con órdenes
negativos — un término descartado por "orden demasiado alto" puede, al multiplicarse después por un
factor de orden NEGATIVO, volver a caer dentro de la ventana de órdenes que sí importan, dando un
resultado mal. **Corregido**: no truncar nunca en `__mul__` (las series base ya vienen truncadas UNA
vez por sympy antes de convertirse a `LSeries`, así que todo lo demás es aritmética EXACTA de
polinomios de Laurent finitos — no hace falta truncar más, y el tamaño no explota: el número de
órdenes distintos en cualquier producto está acotado por la suma de los rangos de los factores).

## Verificación byte a byte contra el pipeline sympy ya confiable (grado9 y grado11, familia 1)

- `fast_grado9_familia1.py`: B1,B2,B3 recalculados con `lseries.py`, comparados coeficiente a
  coeficiente contra `B123.pkl` (ya calculado con sympy, parte 27) — **coinciden exacto**. Tiempo:
  orbit sum 6.5s (vs. varios minutos con sympy).
- `fast_grado11_familia1.py`: S1..S15 recalculados, comparados contra `S115.pkl` — **coinciden
  exacto** en los 15. Tiempo: orbit sum 60s (vs. 321s con sympy, ~5.4x).

## Segundo bug real encontrado (en la familia 2, NO en `lseries.py`): parámetros que dejaban un
## `sqrt(3)` sin cancelar hasta la suma final sobre S₅

Al correr `fast_grado11_familia2.py` contra `S115_familia2.pkl` (parte 32), aparecieron mismatches.
Diagnosticado por descomposición término a término (no aceptado sin explicar): los parámetros de la
familia 2 usados en la parte 32 (`t0=-1, m0sq0=1`) hacían que `R(x)` tuviera un coeficiente líder
irracional (`sqrt(3)`, ya que `sqrt(4/3)=2/√3`) — las polarizaciones INDIVIDUALES (`eps_ser[3],
eps_ser[4],eps_ser[5]`, antes de la suma sobre las 120 permutaciones) genuinamente contienen términos
con `sqrt(3)` que sólo cancelan al simetrizar sobre toda la órbita de S₅ (confirmado: `S115_familia2.pkl`
NO tiene irracionales, el pipeline sympy original los manejó bien todo el tiempo). El bug estaba en
`_to_fraction` (dentro de `lseries.py`): usaba `int()` sobre una expresión sympy sin verificar que
fuera un entero exacto — `int(sqrt(3))` da `1` **silenciosamente**, sin error, corrompiendo el dato
en vez de fallar. **Corregido**: `_to_fraction` ahora exige explícitamente que numerador y
denominador sean enteros sympy, y lanza un error claro si no (en vez de truncar mal sin avisar).

**Corrección de fondo**: se cambiaron los parámetros de la familia 2 a `t0=-3, m0sq0=1` (encontrados
por búsqueda chica, igual que se hizo para la familia 1: "elegir t0, m0sq perfectos
cuadrados/simples para mantener todo en Q") — con esto `sqrt(R(x))` es exactamente racional, sin
ningún término irracional en ningún paso intermedio. Se reconfirmó que el patrón de canales
(9 crecientes / 1 fijo, `(1,3)` el único fijo) es el mismo que antes. Se regeneró `S115_familia2.pkl`
con el pipeline sympy (referencia lenta, ~520s) y se re-verificó `fast_grado11_familia2.py` contra
él: **coinciden exacto en los 15** (orbit sum: 17s vs. 521s, ~30x).

## Resultado de la parte 32, RECONFIRMADO con los datos corregidos

- Los 9 vectores de la familia 1 (parte 31) siguen violando `s²` en la familia 2 (caen en `x⁻¹¹`).
- El subespacio universal (combinaciones que cumplen `≤s²` en AMBAS familias) sigue teniendo
  **dimensión exactamente 2**, con los 2 vectores arrancando exacto en `x⁻⁴=s²` en ambas trayectorias
  — mismo resultado cualitativo que con los parámetros (accidentalmente irracionales, pero
  correctamente manejados por sympy) de la parte 32.

## Verificación de robustez pendiente (marcada como "la más urgente" en la parte 32) — HECHA,
## gracias al backend rápido

Con el pipeline validado, se corrió una **segunda configuración de parámetros independiente en
ambas familias** (`fast_segunda_config.py`, ~2 minutos en vez de ~20-30 minutos que hubiera tomado
con sympy):
- Familia 1: `t0=-4, m0sq=9, cos_phi0=5/13, sin_phi0=12/13` (la misma segunda config ya usada en la
  parte 30 para el piso `s³` de grado 9, reutilizada acá para grado 11).
- Familia 2: `t0=-12, m0sq0=1, cos_phi0=5/13, sin_phi0=12/13` (mismo `m0sq0` que garantiza `sqrt(R)`
  racional, pero `t0` y `phi0` distintos).

**Resultado**: las relaciones exactas `A(S8)=-1/2·A(S1)` y `A(S14)=A(S15)` (parte 30) se reproducen
EXACTAS en esta configuración también. Y, lo más importante: **el subespacio universal vuelve a dar
dimensión exactamente 2, con ambos vectores arrancando en `x⁻⁴=s²` en las dos familias** —
reproducido IDÉNTICO en una segunda configuración de parámetros completamente independiente.

## Estado honesto actualizado

El hallazgo de la parte 32 (subespacio de dimensión 2, grado 11, que satisface `≤s²` en 2 familias
de trayectorias de Regge cualitativamente distintas) ahora está verificado en: (a) 2 familias de
trayectoria distintas (no relacionadas por reetiquetado de S₅), Y (b) 2 configuraciones de
parámetros numéricos independientes dentro de cada familia — el mismo nivel de robustez que el
proyecto exige para los demás resultados centrales (ver partes 27, 30). Sigue pendiente (menor
prioridad, ver parte 32): una tercera familia de trayectoria (p.ej. otro exponente jerárquico, o
multi-Regge con Sudakov) y acotar mejor la dimensión real del espacio de grado 11 (hoy `≥15`).

## Archivos

- `scripts/lseries.py` — backend rápido (Fraction puro, sin sympy).
- `scripts/test_lseries.py` — test de correctitud (100 casos aleatorios vs. sympy).
- `scripts/fast_grado9_familia1.py`, `fast_grado11_familia1.py`, `fast_grado11_familia2.py` —
  versiones rápidas con verificación byte a byte contra los `.pkl` de referencia.
- `scripts/fast_segunda_config.py` — segunda configuración de parámetros en ambas familias
  (robustez), reproduce todo exacto.
- `scripts/regge_n5_familia2_kinematics.py` — corregido a `t0=-3, m0sq0=1` (antes `t0=-1`, que dejaba
  `sqrt(3)` sin cancelar hasta el final).
- `scripts/S115_familia2.pkl`, `scripts/grado11_base_universal.pkl` — regenerados con los parámetros
  corregidos.
