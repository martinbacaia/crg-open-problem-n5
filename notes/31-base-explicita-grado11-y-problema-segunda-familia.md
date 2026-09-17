# Fase 2 — Paso 31: base explícita de los 9 vectores de grado 11 (≤s²) extraída y verificada;
# obstrucción estructural real encontrada al intentar construir una segunda familia de Regge

## Punto 1 (hecho): base explícita de la familia de dimensión 9

Se corrió `scripts/regge_n5_simbolico_grado11.py` (no estaba guardado el `.pkl` de una sesión
anterior) para regenerar `S115.pkl` (las 15 series de Laurent exactas de S1..S15 en la trayectoria
de Regge ya usada, partes 22-23/28/30). Con `scripts/grado11_extraer_base9.py` se construyó la
matriz `6×15` (filas = coeficientes de `x^-10` a `x^-5`) y se extrajo su núcleo exacto (sympy,
aritmética racional): **dimensión 9**, confirmando exactamente la tabla de la parte 30 (rango
6→dim. núcleo 9 en `x^-5`; rango sube a 7 —núcleo 8— al incluir `x^-4`, también coincide).

Se guardaron los 9 vectores explícitos (`grado11_base9.pkl`, cada uno como combinación lineal
racional exacta de `S1..S15`) y se verificó — reconstruyendo cada combinación y buscando el primer
orden no nulo de su serie de Laurent — que **los 9 vectores arrancan exactamente en `x^-4` (`=E⁴=s²`)**,
ninguno antes: confirma de forma independiente (reconstrucción directa, no solo álgebra de rango) que
los 9 satisfacen el bound `≤s²` en esta trayectoria. Los coeficientes son fracciones racionales
grandes (la base que da sympy vía eliminación con pivotes en `S10..S15`, no una base "linda") — útil
para evaluación pero no ilustrativa por sí misma; si se necesita una base más legible habría que
buscar combinaciones enteras pequeñas dentro de este espacio de dimensión 9 (no se intentó, no es
necesario para el punto 3).

## Punto 2 (bloqueado — hallazgo importante antes de seguir)

Se intentó diseñar la segunda familia de límites de Regge pedida (ej. `{1,2}` duro vs. `{3,4,5}`
retroceso). Análisis (sin escribir código, análisis cinemático a mano) reveló una obstrucción real:

**Cualquier reetiquetado de la MISMA construcción de familia 1** (elegir un tríos distinto como
"cluster duro" y el par restante como "retroceso", p.ej. duro=`{1,2,4}`/retroceso=`{3,5}`) da
**resultados numéricamente IDÉNTICOS** a la familia ya evaluada — no es un chequeo nuevo. Razón: cada
`Sᵢ` (y por tanto cualquier combinación lineal, incluidos los 9 vectores del punto 1) se construye
sumando explícitamente sobre las 120 permutaciones de las 5 etiquetas (`S5`-simetrización completa,
ver el código en `regge_n5_simbolico_grado11.py`). Evaluar un objeto `S5`-simétrico en una trayectoria
reetiquetada da, por construcción, el mismo valor que evaluarlo en la trayectoria original. Se
verificó también que construir un límite genuinamente `{1,2}` duro + `{3,4,5}` retroceso (3 cuerpos)
de forma NO degenerada (sin que quede una partícula en reposo o con masa creciente sin límite) requiere,
por conservación de momento, que al menos 2 de las 3 partículas de "retroceso" formen a su vez un
sub-par que crece en energía como boost — es decir, la construcción colapsa estructuralmente a ser
una reetiqueta de la familia `{hard-trío, recoil-par}` ya usada (dim. de "quién es duro" no es un
grado de libertad real distinto, dada la simetría completa de los invariantes que se están probando).

**Conclusión honesta**: una segunda familia que sea una prueba genuinamente independiente necesita
un patrón de "qué crece / qué se mantiene fijo" cualitativamente distinto (no alcanzable por
reetiquetado) — candidatos serios: (a) un límite jerárquico de 2 escalas (ej. dentro del retroceso
`{4,5}` de la familia 1, enviar también su propia energía interna a infinito a una tasa distinta,
`E'~E^a` con `a≠1`, dando un límite de Regge "iterado"/multi-escala genuino); (b) cinemática tipo
multi-Regge (ordenamiento de rapidities) con un patrón de invariantes fijos que no son simplemente los
4 pares "canal fijo" de la familia 1 sino combinaciones de momento transverso — más fiel a la
literatura pero requiere más maquinaria (variables de Sudakov) y no es obviamente el objeto correcto
para testear el bound `s²` en el sentido de CRG (que pide UN invariante creciendo con el resto
literalmente fijo, no un ordenamiento jerárquico de varios). No se implementó ninguna de las dos
todavía — se prefirió reportar la obstrucción y las opciones antes de comprometerse a una
construcción que podría no ser una prueba real.

## Próximo paso

Decisión pendiente del usuario: elegir entre (a) el límite jerárquico de 2 escalas dentro de la
familia ya validada (extensión más simple y directamente ligada a lo ya construido), (b) intentar
multi-Regge con variables de Sudakov (más fiel a la literatura, más trabajo), o (c) otra idea. Ver
mensaje de cierre de esta sesión para el detalle de las opciones.

## Archivos

- `scripts/grado11_extraer_base9.py` — extrae y verifica los 9 vectores.
- `scripts/grado11_base9.pkl` — los 9 vectores (coeficientes exactos en S1..S15).
- `scripts/S115.pkl` — regenerado (series de Laurent exactas de S1..S15, config. de referencia).
