# Fase 2 — Paso 19: hacia el chequeo de Regge growth en n=5 — la definición exacta del paper
# (n=4) y una sutileza real que bloquea una generalización ingenua

## Lo que dice exactamente el paper (releído a fondo, secciones 1.4, 1.5, 2.9-2.10 de
## arXiv:1910.14392)

**La CRG conjecture, tal como se usa en el paper (no la derivación desde chaos bound — eso es el
"Problema A" ya descartado en `fase1-mapa/04-profundizacion-crg.md`)**: cualquier S-matrix físicamente
aceptable de una teoría clásica debe crecer, en el límite de Regge (`s→∞` a `t` FIJO), **no más
rápido que `s²`**. Cualquier interacción de contacto (vértice local, sin polos) que en espacio plano
crecería más rápido que `s²` a `t` fijo se descarta.

**Cómo se aplica concretamente para n=4** (sección 2.10): el módulo de S-matrices polinomiales de 4
fotones está generado por elementos que transforman en las representaciones `1_S`, `1_A`, `2_M` de
S₃ (la misma teoría de representaciones que usamos para validar el caso n=4 al principio de esta
sesión — parte 1). Para cada generador, sus descendientes (multiplicar por potencias de los 2
invariantes S₃-invariantes `stu` y `s²+t²+u²`) tienen una tasa de crecimiento en el límite de Regge
que se puede **derivar exactamente** (ecuaciones 2.59-2.66 del paper) porque se sabe cómo escalan
`stu` y `s²+t²+u²` cuando `s→∞` a `t` fijo (con `u=-s-t→-∞` también, ya que `s+t+u=0`). El resultado
concreto: a orden 2, 4 y 6 derivadas hay estructuras que crecen como `s⁰,s,s²` respectivamente
(saturando el bound justo en 6 derivadas); **todo lo de orden 8 derivadas o más crece
necesariamente más rápido que `s²`** — de ahí que el candidato de 4 parámetros (ecuación 1.6, hasta
6 derivadas) sea exactamente el límite de lo permitido.

**Punto clave para nuestro proyecto**: esta derivación (ecs. 2.59-2.66) es **específica de la
estructura de S₃ actuando en 3 variables (s,t,u)** — no es una fórmula genérica que se pueda aplicar
mecánicamente a n=5. Habría que re-derivarla desde cero para la estructura de S₅ actuando en 5
invariantes de Mandelstam (partes 2-8 de esta sesión), y ni siquiera el paper original lo intenta
para n≥5 (confirmado ya en la Fase 1).

## Antes de eso: hace falta definir qué ES el "límite de Regge" para 5 partículas — no está en la
## literatura, y una generalización ingenua NO funciona

Para n=4, "s→∞ a t fijo" tiene una realización concreta a nivel de VECTORES de momento (no solo de
los símbolos s,t): es el límite de alta energía + ángulo de scattering hacia adelante (forward),
donde LAS 4 partículas participan de un boost correlacionado — no es que 2 partículas se aceleren y
las otras 2 queden quietas.

**Se intentó una generalización ingenua y se encontró que NO funciona, por una razón estructural
real** (no un error de cómputo): la construcción obvia — boostear las partículas 1,2 en direcciones
de cono de luz opuestas (`p₁ ~ Λ e₊`, `p₂ ~ Λ e₋`, con `e₊,e₋` nulos de referencia) mientras se
mantienen fijas 3,4,5 — falla porque:
1. Si p1,p2 crecen, por conservación de momento **algo más también tiene que crecer** (no se puede
   tener `p1+p2` creciendo con `p3+p4+p5` fijo, ya que su suma debe ser exactamente cero) — hace
   falta que al menos una partícula "espectadora" absorba el retroceso, creciendo también con Λ.
2. Peor: **ningún momento nulo genérico puede tener overlap cero con AMBAS direcciones `e₊,e₋`**
   (un vector nulo con componente cero en ambas direcciones de cono de luz viviría enteramente en el
   subespacio transverso, que es espacio-tipo — un vector nulo ahí solo puede ser el vector cero).
   Esto significa que invariantes como `s₁₃` (que en el límite ingenuo uno querría mantener FIJOS,
   análogos a "t") **necesariamance crecen también** si p3 es un momento nulo genérico — contradice
   la idea de "t fijo".

**Conclusión honesta**: en el caso n=4 real, el límite "t fijo" funciona porque las 4 partículas
participan de una estructura muy específica y correlacionada (el boost + forward limit combinados),
no porque 2 partículas se queden quietas. Generalizar esto correctamente a n=5 (dónde "fijar t" debe
significar algo preciso con 5 momentos) es una **decisión de modelado real, no resuelta en la
literatura para n≥5**, y no es evidente cuál es la elección correcta o si hay una sola.

## Por qué no se apuró un chequeo de todos modos

Se decidió **no** construir una cinemática de Regge ad hoc y reportar "pasa" o "no pasa" el bound —
eso arriesgaría publicar una conclusión sobre el objetivo final del proyecto basada en una
definición arbitraria y posiblemente incorrecta del límite. Dado que el objetivo explícito del
proyecto es "publicar algo correcto, no rápido" (y esta sesión ya encontró y corrigió varios
resultados propios erróneos que casi se dan por buenos — partes 2,4,6,14, 18-19), no vale la pena el
riesgo de apurar este paso, que es además el más importante de todo el proyecto (el chequeo final).

## Camino recomendado para la próxima sesión

1. **Primero, replicar el límite de Regge de n=4 a nivel de VECTORES explícitos** (no solo de los
   símbolos s,t,u) — construir la parametrización de cono de luz concreta que realiza "s→∞ a t
   fijo" con los 4 momentos físicos de verdad, y verificar numéricamente que el candidato ya conocido
   del paper (ecuación 1.6, con sus 4 términos) efectivamente satura/respeta el bound `s²` en esa
   realización explícita. Esto sirve de validación (igual que se hizo con la representación de S₃
   antes de atacar S₅) y aclara exactamente qué estructura de cono de luz generalizar.
2. **Después, generalizar esa construcción específica a 5 partículas**, resolviendo con cuidado el
   problema de conservación de momento (necesariamente alguna combinación de partículas "espectadoras"
   absorbe el retroceso) y decidiendo de forma explícita y justificada qué significa "t fijo" con 5
   momentos (candidato natural a explorar: estructura tipo "multi-Regge kinematics", bien establecida
   en la literatura de amplitudes de gauge — cadena de partículas ordenadas por rapidez, con
   invariantes entre vecinos creciendo y el resto acotado).
3. Solo con esa definición en mano, evaluar el candidato de grado 9 ya encontrado (parte 20) y
   determinar su tasa de crecimiento real.

## Estado

Sin cambios en los resultados ya confirmados (grados 7,8,9,10,11 de la tabla de la parte 21). Este
paso es preparatorio: se despejó la definición exacta del bound CRG y se identificó, con una razón
estructural concreta (no solo "parece difícil"), por qué la generalización obvia del límite de Regge
a n=5 no es trivial — se necesita una sesión dedicada a construirlo correctamente antes de aplicar
cualquier chequeo numérico.

## Archivos

- No se generó código nuevo esta parte — el trabajo fue de lectura de fuente primaria
  (`fuentes/classifying-4photon-4graviton-smatrices.pdf`, secciones 1.4, 1.5, 2.9, 2.10) y un intento
  de construcción de cinemática de Regge (hecho a mano/analíticamente, sin llegar a implementarse en
  código, precisamente porque reveló el problema estructural antes de llegar a esa etapa).
