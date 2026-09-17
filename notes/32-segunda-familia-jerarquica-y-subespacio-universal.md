# Fase 2 — Paso 32: segunda familia de límites de Regge (jerárquica, 2 escalas) construida y
# evaluada — los 9 vectores de la parte 30/31 NO sobreviven, pero aparece un subespacio universal
# de dimensión 2 que SÍ satisface `≤s²` en ambas trayectorias

## Motivación y elección de la segunda familia

Tras el hallazgo de la parte 31 (cualquier reetiquetado de la familia 1 da resultados idénticos,
porque `S1..S15` ya están simetrizados sobre las 120 permutaciones de `S₅`), se necesitaba un patrón
de crecimiento/fijeza cualitativamente distinto, no alcanzable por reetiquetado. El usuario eligió
la opción de **límite jerárquico de 2 escalas**: reutilizar exactamente la construcción de la
familia 1 (cluster duro `{1,2,3}`, `s₁₂→∞`, `t₁₃` fijo vía el ángulo ya derivado; retroceso `{4,5}`
boosteado desde un decaimiento de 2 cuerpos) pero, en vez de mantener `M₄₅²` FIJO (como en la
familia 1), se lo hace crecer también con `E`, a la MISMA tasa que `s₁₂` pero como una escala
independiente: `M₄₅² = (m0sq0)²·E²` (en la variable `x=1/E`, `m0sq = m0sq0²/x²`).

**Por qué esto es tratable exactamente**: la derivación algebraica de `cos_theta, E3` (que fija
`t₁₃=t0` y `M₄₅²=m0sq` simultáneamente) es válida para CUALQUIER valor de `m0sq`, incluso dependiente
de `x` — no hubo que rederivar nada, sólo sustituir el valor. Eligiendo el exponente par más simple
(`m0sq0` en vez de `m0sq0²` bajo la raíz, para que `m0=sqrt(m0sq)=m0sq0/x` sea racional en `x`, sin
potencias fraccionarias) se mantiene toda la aritmética `QExt` exacta sin cambios estructurales.

**Detalle técnico encontrado y corregido**: con `m0sq0=2` la fórmula `E3=E-m0sq/(4E)` se anula
idénticamente (`m0sq0²/4=1` cancela exacto contra el `1` de `E`), dando `cos_theta=zoo` (división por
cero) — detectado inmediatamente (todos los checks daban `nan`), diagnosticado (no se aceptó el `nan`
sin explicarlo) y corregido usando `m0sq0=1` (evita la cancelación accidental).

## Verificación de la nueva cinemática

- Masa nula de las 5 partículas: exacto, 0 en toda la serie truncada.
- Conservación de momento: exacto, 0.
- Clasificación de los 10 pares (orden líder en `x`): **sólo `(1,3)` (=`t₁₃`) se mantiene fijo**
  (orden `x⁰`) — los otros 9 pares, incluidos `(2,4),(2,5),(4,5)` (que en la familia 1 eran "canal
  fijo"), ahora **crecen como `x⁻²` (`~E²~s`)**. Esto confirma que es un patrón cualitativamente
  distinto (9 crecientes/1 fijo, vs. 6 crecientes/4 fijos en la familia 1) — no una versión
  reetiquetada de la familia 1.
- La asignación de "compañero creciente" para construir las polarizaciones (`growing_partner`,
  ya validada en la familia 1) sigue siendo segura en la familia 2, porque los pares que usa como
  referencia ya eran crecientes en la familia 1 y **siguen siéndolo** en la familia 2 (verificado
  explícitamente antes de confiar en la construcción, repitiendo la lección de la parte 21/23: nunca
  usar un par de "canal fijo" como referencia de gauge — acá el único fijo es `(1,3)`, que no se usa
  como referencia para ninguna partícula).

## Resultado 1: los 9 vectores de la parte 31 NO sobreviven en la familia 2

Se recalcularon `S1..S15` (grado 11) en la familia 2 (`scripts/regge_n5_familia2_grado11.py`) —
cada `Sᵢ` individual crece como `x⁻¹¹` (`~E^11~s^5.5`, un poco distinto del `~s⁵` de la familia 1,
esperable porque el patrón de canales creciente/fijo es distinto). Se evaluaron los 9 vectores ya
extraídos (`grado11_base9.pkl`) en esta nueva trayectoria: **los 9 caen exactamente en `x⁻¹¹`, sin
ninguna cancelación** — violan el bound `s²` por un margen grande. Es decir, la cancelación que
producía `≤s²` en la familia 1 es una propiedad de ESA trayectoria específica, no universal.

## Resultado 2 (hallazgo positivo): un subespacio de dimensión 2 SÍ satisface `≤s²` en AMBAS familias

Siguiendo la recomendación del propio plan de esta sesión (punto 4: "diagnosticar... podría ser que
sólo un subespacio más chico sea universal"), se construyó el sistema combinado: la matriz de
coeficientes de Laurent de `S1..S15` en AMBAS trayectorias (familia 1: órdenes `x⁻¹⁰..x⁻⁵`, 6 filas;
familia 2: órdenes `x⁻¹¹..x⁻⁵`, 7 filas — 13 filas en total, exigiendo que la combinación se anule en
todos esos órdenes en las 2 trayectorias simultáneamente) y se calculó su núcleo exacto (`sympy`,
racional): **rango 13 (completo), núcleo de dimensión exactamente 2**.

Se verificó DIRECTAMENTE (no sólo por rango/núcleo, sino reconstruyendo cada combinación y buscando
su primer orden no nulo en cada una de las 2 series de Laurent por separado): **los 2 vectores `U1,
U2` arrancan exactamente en `x⁻⁴` (`=s²`) tanto en la familia 1 como en la familia 2** — ninguno
antes, en ninguna de las 2 trayectorias. Guardados en `scripts/grado11_base_universal.pkl`.

## Estado honesto — qué prueba esto y qué NO

**Lo que está probado (exacto, en las 2 configuraciones concretas usadas)**:
- Existe un subespacio de dimensión (al menos) 2 dentro del espacio de grado 11 (dimensión ≥15)
  cuyo crecimiento es `≤s²` simultáneamente en 2 familias de trayectorias de Regge cualitativamente
  distintas (patrones de canal creciente/fijo distintos, no relacionadas por reetiquetado de `S₅`).
- Esto es **la primera evidencia, dentro de este proyecto, de un candidato de grado 11 que sobrevive
  más de una familia de trayectoria** — más débil que el resultado "dimensión 9" de la parte 30 (que
  sólo se había probado en 1 familia), pero genuinamente más sólido en el sentido que importa (no
  depende de una elección particular de qué se mantiene fijo).

**Lo que NO está probado todavía**:
1. **Una sola configuración de parámetros por familia**: no se repitió (todavía) con una segunda
   configuración numérica de `t0,m0sq0,phi0` dentro de cada familia — el propio patrón de esta
   sesión (y de sesiones anteriores) es que un resultado "sobrevive 1 configuración" puede no
   sobrevivir una segunda (ver partes 27/30, donde SIEMPRE se pidió una segunda config antes de
   confiar). **Este es el paso más urgente pendiente antes de anunciar nada.**
2. Sólo se probaron 2 familias (de las posiblemente muchas generalizaciones razonables de "límite de
   Regge para n=5"). El bound CRG, para sostenerse genuinamente, debería sobrevivir CUALQUIER
   generalización razonable — 2 no es una muestra grande, aunque son cualitativamente distintas
   (patrón de canales 6/4 vs. 9/1).
3. No se determinó si `U1,U2` (o alguna combinación de ambos) puede bajar por debajo de `s²` en
   alguna de las 2 trayectorias (no se buscó, no es necesario para el punto 3 del plan de la sesión,
   pero sería información extra).
4. La dimensión real del espacio de grado 11 podría ser mayor a 15 (parte 28) — si aparecieran más
   generadores independientes, el subespacio universal de dimensión 2 podría no ser exhaustivo (podría
   haber más direcciones universales fuera de las 15 ya conocidas).

## Próximo paso concreto

1. **Prioridad más alta**: repetir la evaluación de `U1,U2` con una segunda configuración de
   parámetros en CADA familia (análogo a como se hizo en la parte 30 para la familia 1 sola) — si el
   resultado `≤s²` en ambas familias se sostiene también ahí, sería mucho más sólido.
2. Extraer una forma más legible de `U1,U2` (los coeficientes actuales son fracciones grandes de la
   eliminación de sympy) — buscar combinaciones enteras pequeñas dentro del mismo subespacio de
   dimensión 2.
3. Considerar una TERCERA familia (p.ej. variando el exponente jerárquico, `m0sq~E^b` para otro `b`,
   o explorando multi-Regge con variables de Sudakov) para seguir acotando qué tan universal es este
   subespacio de dimensión 2.

## Archivos

- `scripts/regge_n5_familia2_kinematics.py` — cinemática de la familia 2 (jerárquica, 2 escalas).
- `scripts/regge_n5_familia2_eps.py` — polarizaciones (misma referencia de gauge que familia 1,
  reverificada segura).
- `scripts/regge_n5_familia2_classify.py` — clasificación de los 10 pares (9 crecientes/1 fijo).
- `scripts/regge_n5_familia2_grado11.py` — `S1..S15` en la familia 2 (→ `S115_familia2.pkl`).
- `scripts/grado11_evaluar_base9_familia2.py` — evalúa los 9 vectores de la parte 31 en familia 2
  (resultado: los 9 violan `s²`).
- El cálculo del subespacio universal (`grado11_base_universal.pkl`, dimensión 2) se hizo con un
  script ad hoc de pocas líneas (sympy `Matrix.nullspace()` combinando ambas familias) — reproducible
  en <1 minuto a partir de `S115.pkl` y `S115_familia2.pkl`.
