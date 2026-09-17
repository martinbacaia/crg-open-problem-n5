# Fase 2 — Paso 28: grado 11 recalculado con sandwich — dimensión ≥15 (no 3 como en la parte 21)

## Resultado

Siguiendo la advertencia de consistencia dejada en la parte 27 ("si grado 9 estaba incompleto por no
incluir el sandwich, grado 11 probablemente también"), se recalculó grado 11
(`scripts/grado11_completo_con_sandwich.py`) agregando sistemáticamente una nueva familia de
semillas:

- Familia A (ya usada en la parte 21): 165 semillas.
- Familia B (ya usada en la parte 21, no exhaustiva): 2 semillas.
- **Familia C (NUEVA)**: "1 pareja F:F + 1 sandwich + 1 singlete C" × corrector Mandelstam grado 2
  (grado 2+4+3+2=11) — 6 formas de contraer el sandwich × 6 formas de contraer el singlete × 4
  correctores representativos = 144 semillas.

**Total: 311 semillas (familias A+B+C)**, evaluadas con `mpmath` a 60 dígitos (lección de la parte
21: float64 no es confiable en grado 11) en 2 lotes independientes de 15 y 20 muestras cinemáticas.

**Resultado: rango = 12 en ambos lotes**, con residual tras el último pivote de `1.7-2.0×10⁻⁶⁰` —
esencialmente exacto.

**Se agregó además la familia D (NUEVA)**: "2 sandwiches + 1 singlete C" (grado 4+4+3=11, SIN
corrector — topología sin análogo en grado 9, ahí no entraban 2 sandwiches) — 27 semillas
representativas (`scripts/grado11_con_familia_D.py`). **Con las 4 familias (A+B+C+D, 338 semillas
totales), el rango sube a 15** — verificado estable en 2 lotes independientes (20 y 25 muestras,
seeds 1-20 y 201-225), con residuales `2.9×10⁻⁶¹` y `2.5×10⁻⁶¹` — esencialmente exactos.

**Esto NO es necesariamente la dimensión final** — las familias C y D usadas son submuestras
razonables, no exhaustivas (se restringieron a subconjuntos representativos de correctores y
contracciones en vez de la combinatoria completa). La dimensión REAL de grado 11 es **al menos 15**,
posiblemente más.

## Por qué esto es tan distinto de la parte 21 ("dimensión 3")

La parte 21 solo probó topologías basadas en `F^i:F^j` y `C^k_{ab}` — nunca incluyó el sandwich
(igual que grado 9 en la parte 20). El sandwich añade una familia completamente nueva de
contracciones (`M^{ij}_{ab}` con múltiples elecciones de índices) que resulta en MUCHAS más
direcciones independientes de las que se sospechaba. Esto es consistente con lo que ya se vio en
grado 9 (de dimensión "1" a "3"), pero la escala del salto en grado 11 (de 3 a ≥12) es bastante
mayor — sugiere que la topología sandwich se vuelve combinatoriamente más rica a medida que aumenta
el grado (tiene más formas de contraerse con más momentos disponibles).

## Consecuencia inmediata

Con dimensión ≥15 (mucho mayor que se pensaba), el margen para buscar combinaciones que satisfagan
Regge growth es **mucho más grande** de lo que parecía en la parte 24 (donde con dimensión 3 se
concluyó que el piso era `s⁴`, igual que grado 9). **Esa conclusión de la parte 24 queda descartada**
— no se puede confiar en un análisis de cancelación hecho sobre una familia incompleta (dimensión 3
de 15+).

## Qué falta (próximo paso concreto)

1. Repetir el análisis de cancelación en el límite de Regge (parte 23) con la base completa de 15
   direcciones independientes ya encontradas — dado el precedente de grado 9 (2 de 3 direcciones
   bastaron para bajar de `s⁴` a `s³`), con 15 direcciones en grado 11 hay bastante más espacio para
   buscar una combinación que baje significativamente, quizás hasta `s²`.
2. Evaluar si vale la pena seguir ampliando la búsqueda de generadores en grado 11 (la combinatoria
   completa de las familias C y D, no solo los subconjuntos representativos ya probados) antes o
   después del paso 1 — dado que ya hay bastante margen con 15, puede ser más eficiente ir directo
   al análisis de Regge growth primero.

## Archivos

- `scripts/grado11_completo_con_sandwich.py` — familias A+B+C, rango 12 confirmado en 2 lotes.
- `scripts/grado11_con_familia_D.py` — familias A+B+C+D, rango 15 confirmado en 2 lotes
  independientes adicionales.
