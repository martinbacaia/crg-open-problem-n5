# Fase 2 — Paso 37: la topología nueva T ("sandwich triple") NO rescata el bound de Regge en
# grado 11 — el núcleo combinado de las 3 familias sigue siendo `{0}`, ahora probado en el espacio
# COMPLETO de dimensión 16 (no sólo en las 15 direcciones viejas)

## Motivación

Pregunta pendiente al cierre de la parte 36: la topología T (sandwich triple, `T(i,j,k;a,b) :=
p_a.F^i.F^j.F^k.p_b`) sube la dimensión real de grado 11 de 15 a 16 (dentro de las topologías
exploradas) — pero nunca se había evaluado esa dirección nueva en el límite de Regge. La pregunta
concreta (marcada como el paso de mayor valor pendiente): ¿la dirección nueva de T, sola o
combinada con las 15 anteriores, permite que el subespacio universal (que satisface `≤s²` en las 3
familias de Regge ya construidas — partes 22-23, 32, 34) deje de ser `{0}`?

## Implementación

1. **`Tsandwich` verificado en el backend simbólico** (`tsandwich_cached_check.py`): se derivó una
   versión de `T` basada en las mismas cachés `eta_pp/eta_pe/eta_ee` que ya usa `Msandwich` en el
   pipeline de Regge (en vez de la cadena genérica sobre vectores explícitos de
   `topologia_nueva_triple_sandwich.py`, más lenta) — verificada EXACTA contra la referencia en 3
   semillas cinemáticas (2820 tuplas `(i,j,k,a,b)` cada una, coincidencia perfecta, `mpmath` 60
   dígitos). Necesario antes de confiar en `T` dentro del pipeline de Laurent.
2. **`S16` definido**: la variante más limpia de T+C+C encontrada en la parte 36 —
   `S16(i,j,l,m,k) := T(i,j,k;l,m)·C(l;i,j)·C(m;i,k)`, simetrizada sobre las 120 permutaciones de
   S₅ (mismo patrón que S1..S15). Implementado en `basis_S116.py` (backend rápido `lseries.py`,
   reusable para las 3 familias).
3. **Cross-check obligatorio**: para las familias 1 y 2, S1..S15 recalculados con `basis_S116.py`
   coinciden EXACTO, coeficiente a coeficiente, con `S115_fast.pkl`/`S115_familia2_fast.pkl` (ya
   validados en sesiones anteriores contra el pipeline sympy original) — confirma que el
   refactor a un módulo compartido no introdujo ningún error antes de confiar en `S16`.

## Resultado 1: S16 SÍ agranda el subespacio ≤s² en familias 1 y 2 individualmente

- Familia 1 sola: el núcleo (`≤s²`, órdenes `x⁻¹⁰..x⁻⁵`) crece de dimensión 9 (S1..S15) a
  **dimensión 10** (S1..S16) — S16 aporta una dirección genuinamente nueva ahí.
- Familias 1+2 combinadas (mismo criterio de las partes 32-33): el núcleo crece de **dimensión 2**
  (el subespacio universal ya conocido, con coeficiente de S16 = 0 en sus 2 generadores) a
  **dimensión 3** — el generador nuevo (`U3`) tiene coeficiente de S16 exactamente 1, confirmando
  que es una dirección genuinamente nueva y no una redundancia de las 15 anteriores.

## Resultado 2 (el que responde la pregunta real): la familia 3 mata `U3` igual que mató a `U1,U2`

Se evaluó `U3` en la familia 3 (`M45²~E⁴`, la que ya había matado el subespacio de dimensión 2 en
la parte 34), con el mismo rigor exigido en esa parte (2 truncaciones independientes, `N=30` y
`N=45`, verificadas coeficiente a coeficiente idénticas en todo el rango `x⁻²²..x⁻⁵` usado, no sólo
en el orden líder):

- `U1`, `U2` (las 2 direcciones viejas): arrancan en `x⁻²²` en familia 3 — igual que antes.
- **`U3` (la dirección que involucra a S16): también arranca en `x⁻²²`** — sin ninguna mejora
  respecto de las direcciones viejas.

## Resultado 3 (definitivo): núcleo combinado de las 3 familias, en el espacio COMPLETO de
## dimensión 16, calculado directamente — es `{0}`

En vez de conformarse con evaluar `U1,U2,U3` por separado (que deja abierta la posibilidad de que
alguna combinación de las 3 sobreviva, como advertía la parte 34), se armó la matriz combinada de
restricciones de las 3 familias a la vez sobre las 16 columnas (`S1..S16`): familia 1 (6 filas,
órdenes `x⁻¹⁰..x⁻⁵`), familia 2 (7 filas, `x⁻¹¹..x⁻⁵`), familia 3 (18 filas, `x⁻²²..x⁻⁵`) — 31
filas en total, aritmética racional exacta (`sympy.Rational`, sin ningún redondeo).

**Resultado: rango = 16 (máximo posible), núcleo = `{0}` exacto.** Esto prueba, no sólo sugiere,
que **ninguna combinación lineal de S1..S16 satisface el bound `≤s²` en las 3 familias de Regge
simultáneamente** — la topología T no rescata el candidato de grado 11, ni sola ni en ninguna
combinación con las 15 direcciones ya conocidas.

## Conclusión honesta

La pregunta concreta que quedaba pendiente desde la parte 36 (y era el paso de mayor valor del
proyecto) tiene una respuesta clara y negativa, verificada con el mismo nivel de rigor que el resto
del proyecto (aritmética exacta, 2 truncaciones independientes estables, cross-check contra el
pipeline ya validado). El escenario ya anticipado en la parte 34 ("una tercera familia con un
patrón más agresivo probablemente mata cualquier subespacio pequeño") se confirma también para el
espacio ampliado por T.

**Esto NO cierra la línea de trabajo del grado 11**, pero sí cierra definitivamente esta rama
concreta (la topología T, tal como se construyó). Las opciones que quedan, en orden de expectativa
(igual que se dejó planteado en la parte 36):

1. Buscar un átomo genuinamente distinto de T (no una ampliación de la misma idea) — Levi-Civita
   es el candidato más obvio no explorado todavía.
2. Reconsiderar si grado 11 es el grado correcto — subir a grado 13 (el siguiente grado impar,
   recordar el teorema de paridad de la parte 21) y repetir el programa completo.
3. Repensar si el límite de Regge usado (partes 22-23, 32, 34 — una elección de modelado propia,
   no de la literatura) es la generalización correcta a n=5 de la definición original de CRG
   (limitada a n=4 en el paper de 2019) — dado que 3 familias cualitativamente distintas de esa
   misma elección de modelado ya mataron 2 subespacios candidatos distintos, vale la pena
   cuestionar la elección de modelado misma, no sólo seguir ampliando el espacio de invariantes.

## Archivos

- `scripts/tsandwich_cached_check.py` — verifica la versión cacheada de `T` contra la referencia.
- `scripts/basis_S116.py` — S1..S16 en el backend rápido (`lseries.py`), reusable por las 3
  familias.
- `scripts/fast_grado11_familia1_S16.py`, `fast_grado11_familia2_S16.py`,
  `fast_grado11_familia3_S16.py` — S1..S16 por familia, con cross-check contra S1..S15 ya
  validados.
- `scripts/S116_familia1.pkl`, `S116_familia2.pkl`, `S116_familia3_N30.pkl`,
  `S116_familia3_N45.pkl` — resultados guardados (N=30 y N=45 de familia 3 coinciden exacto en
  todo el rango usado, `x⁻²²..x⁻⁵`).
- `scripts/grado11_S16_nullspace_familia1.py`, `grado11_S16_nullspace_combinado.py` — núcleos
  parciales (familia 1 sola; familias 1+2).
- `scripts/evaluar_U3_familia3.py` — evalúa U1,U2,U3 (núcleo de familias 1+2) en familia 3.
- `scripts/grado11_S16_nullspace_3familias.py` — el cálculo definitivo: núcleo combinado de las 3
  familias a la vez, dimensión 16 → rango 16 → núcleo `{0}`.
