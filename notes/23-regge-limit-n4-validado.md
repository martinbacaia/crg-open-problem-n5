# Fase 2 — Paso 21: límite de Regge de n=4 construido y validado explícitamente a nivel de
# vectores — con un bug real encontrado y corregido que anticipa la dificultad de n=5

## Construcción

Se implementó la realización física estándar del límite de Regge/alta energía a momentum transfer
fijo para scattering elástico 2→2 (`scripts/regge_limit_n4_validacion.py`), en D=6 (consistente con
el resto del proyecto):

```
p1 = -E(1, n1)      p2 = -E(1, -n1)      (par "incoming", convención all-outgoing: negados)
p3 =  E(1, n3)      p4 =  E(1, -n3)      (par "outgoing")
```

con `n1, n3` vectores espaciales unitarios en un plano de scattering, ángulo `θ` entre ellos. Esto da
`s = 4E²` (crece con `E`) y `t = -2E²(1-cosθ)`, que se mantiene **exactamente fijo** en un valor
`t₀` eligiendo `cosθ(E) = 1 + t₀/(2E²)` — es decir, `θ→0` (límite hacia adelante) correlacionado con
`E→∞`, la realización física estándar de "s→∞ a t fijo". A diferencia del intento ingenuo de la
parte 22 (2 partículas boosteadas, 2 quietas — que fallaba), acá **las 4 partículas participan** del
límite de alta energía; lo que se mantiene fijo es el ÁNGULO relativo, no que alguna partícula se
"quede quieta".

## Bug encontrado y corregido: la elección de referencia de gauge no puede usar un socio del
## "canal t"

Al construir explícitamente las polarizaciones `εᵢ` (proyectando un vector semilla fijo para
imponer `εᵢ.pᵢ=0`, técnica usada toda la sesión), el primer intento usó la misma convención de
"otro" que se venía usando en las partes anteriores (una elección más o menos arbitraria de socio).
**Resultado inesperado**: `(F¹:F²)(F³:F⁴)` creció como `Λ⁵`, no `Λ⁴` como predice el conteo
dimensional ingenuo (grado 4 en momentos → crecimiento como `Λ⁴~s²`).

**Diagnóstico** (no se aceptó el resultado raro sin explicarlo): los pares `(1,3)` y `(2,4)` son
precisamente los pares del **"canal t"** — `η(p₁,p₃) = η(p₂,p₄) = -t/2`, que por construcción se
mantiene **fijo** (no crece con `E`), ya que es literalmente el momentum transfer que se está
manteniendo constante. Usar `p₂` como referencia para fijar el gauge de `ε₄` (como hacía la
convención heredada de sesiones anteriores) divide por una cantidad que NO crece, mientras el
numerador sí crece con `E` — la parte transversal de `ε₄` explota espuriamente, sin ningún
significado físico (es un artefacto de la elección de gauge, no del límite físico).

**Corrección**: usar como referencia de gauge el socio del **canal s** (pares `(1,2)` y `(3,4)`, con
`η(p₁,p₂)=η(p₃,p₄)=-s/2`, que sí crece como `E²`) — evita cualquier denominador que se mantenga
chico en el límite. Tras la corrección: `F¹:F² ~ Λ².⁰⁰⁰⁰ ~ s^1.0000` y
`(F¹:F²)(F³:F⁴) ~ Λ⁴.⁰⁰⁰⁰ ~ s^2.0000` — **coincide exacto** con el conteo dimensional ingenuo
esperado (y con la tabla de la sección 2.10 del paper: estructuras genéricas de 2 derivadas crecen
como `s`, de 4 derivadas como `s²`).

## Por qué este bug importa para n=5 (lección directa, no solo curiosidad de n=4)

Este es exactamente el tipo de sutileza que ya se había anticipado (de forma más vaga, "hay que
tener cuidado") en `22-hacia-regge-growth-n5.md`: **en cualquier límite de Regge, hay pares de
partículas cuyo invariante se mantiene deliberadamente fijo (chico) — y usar esos pares como
referencia de gauge para construir polarizaciones explícitas produce artefactos espurios de
crecimiento que no tienen nada que ver con la física real del S-matrix**. Para n=5, con más
invariantes y más estructura de "qué se mantiene fijo vs qué crece" (el problema de conservación de
momento y espectadores ya identificado en la parte 22), este tipo de error será más fácil de cometer,
no más difícil — hace falta identificar con cuidado, ANTES de construir cualquier polarización
explícita, cuáles pares de partículas son "canal creciente" y cuáles son "canal fijo" en la
cinemática de Regge elegida.

## Estado

**Límite de Regge de n=4 validado a nivel de vectores explícitos** — primera vez en el proyecto que
se construye una cinemática de Regge genuina (no solo se manipulan los símbolos s,t,u
simbólicamente). Sirve de plantilla concreta y ya depurada (con su bug ya encontrado y corregido)
para atacar la generalización a n=5 con mayor confianza.

## Qué falta (próximo paso concreto)

1. Generalizar esta construcción (par que crece + ángulo correlacionado) a 5 partículas, resolviendo
   con cuidado el problema de conservación de momento identificado en la parte 22 (con 5 partículas,
   necesariamente alguna(s) deben absorber el retroceso si 2 partículas se separan a alta energía) —
   candidato natural: generalizar a una estructura de 2 "clusters" (p.ej. `{1,2}` de alta energía
   contra `{3,4,5}`, con la subestructura interna de cada cluster fija) y verificar que TODOS los
   invariantes fuera del "canal creciente" elegido se mantengan efectivamente fijos.
2. Identificar explícitamente, en esa nueva cinemática, qué pares son "canal fijo" (para NO usarlos
   como referencia de gauge al construir las εᵢ, repitiendo el error ya encontrado y corregido acá).
3. Evaluar la tasa de crecimiento del candidato de grado 9 ya confirmado (parte 20) en esa
   cinemática, y compararla contra el bound `s²`.

## Archivos

- `scripts/regge_limit_n4_validacion.py` — construcción explícita, con el bug de referencia de
  gauge documentado inline y corregido, y el ajuste log-log final que confirma `s^1` y `s^2`.
