# Fase 2 — Paso 34: tercera familia de Regge (jerarquía más rápida, `M45²~E⁴`) construida y
# verificada — el subespacio universal de dimensión 2 (partes 32-33) NO sobrevive

## Motivación

Con el backend rápido ya validado (parte 33), se atacó el punto pendiente de mayor prioridad: una
TERCERA familia de trayectorias de Regge, genuinamente distinta de las 2 anteriores, para seguir
poniendo a prueba el subespacio universal de dimensión 2.

## Intento fallido documentado (antes de construir nada): variar `t13` directamente rompe la paridad

Primer intento: en vez de escalar `M45²` (como en la familia 2), escalar `t13` mismo a una tasa
intermedia (`t13~E¹`, ni fijo ni `~E²`). **Detectado ANTES de construir nada más** (chequeo de
racionalidad de `sqrt(R)`, la misma disciplina de las familias 1 y 2): esto rompe la paridad de
`R(x)` y da `sqrt(R)` con potencias SEMI-ENTERAS de `x` (`x^(1/2), x^(3/2), ...`), incompatible con
toda la maquinaria (`qext.py`, `lseries.py`) que asume potencias enteras. Diagnosticado: el mecanismo
de "ángulo tunado" que fija `t13` sólo admite escalas de `M45²` con exponente PAR — consistente con
el teorema de paridad ya probado en la parte 21 del proyecto. Se abandonó esta vía sin gastar cómputo
en ella.

## Familia 3 construida: `M45² ~ E⁴` (crece MÁS RÁPIDO que `s₁₂`, no más lento)

Misma estructura que familias 1/2 (`t13` fijo vía el ángulo ya derivado), pero con `M45²` escalando
un exponente par más alto (`m0sq0⁴/x⁴`, es decir `M45²~E⁴`, dos órdenes por encima de `s₁₂~E²`).
Parámetros (`t0=1, m0sq0=1`) encontrados por búsqueda chica para `sqrt(R)` racional (mismo criterio
que familias 1/2). Verificado exacto (álgebra `QExt`, sin truncar): masa nula de las 5 partículas,
`M45²=m0sq` exacto, conservación de momento exacta.

**Clasificación de canales — patrón cualitativamente nuevo**: `(1,3)` sigue fijo (`x⁰`, es `t13`);
`(1,2),(1,4),(1,5)` crecen `~x⁻²` (igual que antes); pero `(2,3),(2,4),(2,5),(3,4),(3,5),(4,5)`
ahora crecen `~x⁻⁴` — **más rápido que el propio canal duro `s₁₂`**. Ningún par de esta lista nueva
es el usado como referencia de gauge (`growing_partner`, reutilizado sin cambios de las familias
1/2) — reverificado explícitamente antes de confiar en las polarizaciones.

## Desafío numérico nuevo (encontrado y resuelto con cuidado, no ignorado)

Como algunos pares crecen mucho más rápido, los invariantes de grado 11 (`S1..S15`) tienen orden
líder mucho más negativo que en las familias 1/2 (`x⁻¹⁹` a `x⁻²²`, contra `x⁻¹¹`). Un primer chequeo
ingenuo (comparar TODOS los coeficientes entre `N=16` y `N=20`, desde el orden líder hasta `x⁰`)
mostró **inestabilidad severa ya en órdenes "subliderantes"** (p.ej. `x⁻⁷`, muy por debajo del orden
líder real `~x⁻²²`) — diagnosticado como el mecanismo de "encogimiento de validez" ya anticipado en
la parte 30 (multiplicar por un factor de orden muy negativo consume más margen de validez de lo
ingenuamente esperado, y acá los factores individuales llegan a `x⁻⁴`, no `x⁻¹` como antes). **No se
aceptó el primer intento sin diagnosticar**: se verificó específicamente que el ORDEN LÍDER (lo único
que hace falta para la pregunta de esta parte) SÍ es estable — coincide exacto entre `N=16` y `N=20`
en los 15 `S_i`, y también al evaluar la combinación `U1,U2` (el subespacio universal de la parte 33)
en esta familia.

## Resultado: el subespacio de dimensión 2 NO sobrevive la familia 3

`U1` y `U2` (universales para familias 1+2, partes 32-33) **arrancan en `x⁻²²` en la familia 3** —
el mismo orden "genérico" que cualquier `S_i` individual, sin ninguna cancelación especial. **Violan
el bound `s²` (que exige orden `≥ x⁻⁴`) por un margen enorme.**

Existe una única combinación (`a·U1+b·U2`, con `a,b` calculados a partir del cociente de los
coeficientes líderes, ambos no nulos) que cancela ese primer orden `x⁻²²` — pero esa dirección
tendría que sobrevivir además TODOS los órdenes subsiguientes de la familia 3 hasta `x⁻⁴` (unos 18
órdenes más de restricciones, cada uno generalmente independiente) usando exactamente 0 grados de
libertad extra (ya se gastó el único que quedaba en las familias 1+2). **No hay ninguna razón
estructural para esperar que sobreviva, y verificarlo exhaustivamente requeriría un `N` de
truncación mucho mayor (con costo computacional creciente más rápido que lineal, ya observado al
intentar `N=30`)** — se decidió NO perseguir ese cálculo dado lo improbable del resultado y el costo
creciente, en vez de gastar cómputo en una dirección de valor esperado casi nulo.

## Conclusión honesta

El subespacio de dimensión 2 encontrado en las partes 32-33 — que sobrevivía 2 familias de
trayectoria (jerárquica `M45²` fijo/`~E²`) y 2 configuraciones de parámetros cada una — **NO es
universal**: una tercera familia, con un patrón de canales genuinamente más agresivo (`M45²~E⁴`), lo
mata de inmediato. Esto es consistente con la expectativa genérica (un espacio de sólo 2 dimensiones
libres difícilmente sobrevive una restricción adicional sustancial) y con el patrón ya visto en este
proyecto de que los resultados "positivos" tienden a ser más frágiles de lo que parecen al principio.
**No se encontró, hasta ahora, ningún candidato de grado 11 (fotones, n=5) que satisfaga el bound
`s²` en más de 2 familias de trayectorias de Regge simultáneamente.**

## Qué queda abierto

1. Explorar si existe ALGUNA combinación de grado 11 universal a las 3 familias, calculando (con un
   `N` de truncación mayor, aceptando el costo) el rango/núcleo combinado de las 3 familias a la vez
   — dado el resultado de esta parte, la expectativa razonable es que ese subespacio sea `{0}`.
2. Reconsiderar si grado 11 es genuinamente el candidato correcto, o si hay que ampliar la búsqueda
   de generadores (dimensión real del espacio hoy `≥15`, parte 28) antes de descartar el grado.
3. Alternativamente, subir a grado 13 (el siguiente grado impar) y repetir el programa completo.

## Archivos

- `scripts/regge_n5_familia3_kinematics.py`, `regge_n5_familia3_eps.py`,
  `regge_n5_familia3_classify.py` — construcción y clasificación de canales de la familia 3.
- `scripts/fast_grado11_familia3.py` — S1..S15 en la familia 3 (backend rápido).
- `scripts/S115_familia3_N16.pkl`, `S115_familia3_N20.pkl` — 2 órdenes de truncación independientes,
  coinciden exacto en el orden líder de los 15 invariantes y de `U1,U2`.
