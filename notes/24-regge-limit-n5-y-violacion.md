# Fase 2 — Paso 23: límite de Regge generalizado a n=5, construido y validado — el candidato de
# grado 9 VIOLA el bound CRG (crece como s⁴, no ≤ s²)

## Construcción del límite de Regge para n=5

Generalizando la construcción validada para n=4 (parte 21), se construyó una familia de cinemáticas
de 5 partículas parametrizada por `E→∞` (`scripts/regge_limit_n5.py`):

- **Sistema "duro"**: `p1, p2` (como en n=4, `s₁₂=4E²→∞`) y una tercera partícula dura `p3`, con
  energía `E3 = E - m₀²/(4E)` y ángulo `θ(E)` respecto de `p1` elegidos para mantener **fijos
  simultáneamente** `t₁₃:=-(p1+p3)²=t₀` y `M₄₅²:=-(p1+p2+p3)²=m₀²` (masa invariante del sistema de
  retroceso).
- **Sistema de retroceso**: `P₄₅:=-(p1+p2+p3)`, un vector tipo tiempo de masa invariante FIJA
  (`m₀`, no crece con `E`) que se aleja con energía creciente — se bombea (boost relativista
  estándar) una configuración de decaimiento de 2 cuerpos hacia `p4, p5`.

**Ajuste encontrado empíricamente para tener crecimiento limpio**: la dirección del decaimiento
`{4,5}` en el frame de reposo de `P₄₅` no puede ser aleatoria genérica (da crecimiento oscilante, sin
ley de potencias definida en varios pares) ni exactamente alineada con el eje del boost (degenera —
varios invariantes se anulan idénticamente, configuración de medida cero). La elección que da
crecimiento limpio en los 10 pares, sin degenerar, es una **inclinación fija pequeña** (`φ₀`)
respecto del eje del boost.

**Clasificación resultante de los 10 pares** (verificada, no asumida, con ajuste log-log exacto):
canal creciente `(1,2),(1,4),(1,5),(2,3),(3,4),(3,5)` — todos `~E²` exactos; canal fijo
`(1,3),(2,4),(2,5),(4,5)` — todos constantes exactas. Se repite la lección de la parte 21: las
polarizaciones εᵢ se construyeron usando SOLO pares del canal creciente como referencia de gauge.

## Resultado: el candidato de grado 9 (parte 20) crece como `s⁴`, violando el bound `s²`

Se evaluó el único invariante S₅-invariante de grado 9 ya confirmado
(`F^i:F^j · F^l:F^m · C^k_{ij} · (d_ik-d_jk)`, simetrizado sobre S₅) a lo largo de esta familia, para
`E` de 50 a 100 000:

```
S_grado9 ~ E^(7.9998)  ~ s12^(3.9999)
```

**Ajuste extremadamente limpio** (pendiente log-log a 4 cifras de precisión, sin dispersión) —
verificado **robusto** en 4 configuraciones independientes de parámetros fijos
(`t₀, m₀², φ₀`) y semillas de polarización distintas, todas dando el mismo exponente
`E^8.00 ± 0.01` (`s^4.00 ± 0.003`):

| `t₀` | `m₀²` | `φ₀` | pendiente en `s₁₂` |
|------|-------|------|---------------------|
| -1.0 | 3.0   | 0.3  | 3.9998 |
| -2.5 | 1.0   | 0.7  | 3.9985 |
| -0.3 | 10.0  | 1.1  | 3.9998 |
| -5.0 | 0.5   | 0.5  | 4.0027 |

**Conclusión: el candidato de grado 9 encontrado en la parte 18-20 (el único invariante
multilineal S₅-invariante, gauge-invariante y Lorentz-invariante de grado mínimo, más allá de las
topologías descartadas de grado 7) crece como `s⁴` en el límite de Regge — excede el bound
conjeturado `≤s²` por 2 potencias, y por lo tanto NO puede ser, por sí solo, parte de una S-matrix
física de fotones en n=5 (si la CRG conjecture generaliza a n=5 de la forma esperada — mismo bound
`s²`, en cualquier canal de 2 partículas).**

## Por qué esto tiene sentido (no es una sorpresa alarmante) y qué significa para seguir

Por conteo dimensional ingenuo, un invariante homogéneo de grado `d` en momentos evaluado en un
límite donde SOLO 2 partículas definen el canal creciente crecería como `Λ^d ~ s^(d/2)`. Para grado
9 (impar), esto da naively `s^4.5` (no entero) — pero el resultado medido es exactamente `s^4`, una
potencia MENOS de lo naive. Esto es coherente: el candidato es una suma sobre 120 términos (órbita de
S₅), y no todos los términos usan pares del "canal creciente" en todos sus factores — el término
dominante en el límite pierde exactamente 1 potencia de crecimiento respecto del máximo teórico
(`E⁹→E⁸`), lo cual convenientemente da una potencia entera limpia de `s`. Es una buena señal de
consistencia interna, no un error.

Esto es análogo a lo que encuentra el paper original para n=4 (sección 2.10): recién los invariantes
de grado 2,4 (y una combinación especial de grado 6) satisfacen el bound `s²`; **todo lo de grado 8 o
más viola necesariamente el bound**. Nuestro grado 9 (bastante más alto que el techo de n=4) fallando
el bound es consistente con ese patrón — no hay razón a priori para esperar que el primer invariante
que aparece en la clasificación (grado 9, bastante alto ya de por sí) sea automáticamente compatible
con Regge growth.

**Esto NO cierra el proyecto ni dice "no hay S-matrix de fotones para n=5"** — dice que la búsqueda
tiene que seguir: (a) hacia invariantes de grado MÁS BAJO que 9 si existieran con otras topologías no
exploradas (poco probable dado el teorema de la parte 19, pero no descartado con el 100% de
certeza), o (b) más realista, hacia **combinaciones lineales de varios invariantes en distintos
grados** cuyos términos de crecimiento máximo se CANCELEN entre sí (exactamente la lógica del caso
n=4: el S-matrix físico de 4 fotones es una combinación de estructuras de grado 2,4,6 con
coeficientes ajustados para que el crecimiento neto no exceda `s²`, no cualquier invariante suelto).
Como el espacio de invariantes en grado 9 es de dimensión 1 (parte 20), no hay margen de cancelación
DENTRO de ese grado — haría falta construir la representación en grados más altos (11, 13, ...) y ver
si combinaciones cruzadas entre esos grados y el de grado 9 (multiplicados por invariantes de
Mandelstam del anillo ya calculado) pueden cancelar el término `s⁴` dominante.

## Advertencia honesta sobre el alcance de esta conclusión

1. Esta es **una elección específica** (aunque bien motivada y validada contra n=4) de generalización
   del límite de Regge a 5 partículas — la literatura no define esto para n≥5. Es razonable, y se
   validó exhaustivamente que reproduce el comportamiento correcto y esperado para n=4, pero no es
   la única generalización posible en principio.
2. Se asumió que el bound CRG generalizado a n=5 sigue siendo `≤s²` (mismo exponente que n=4) en
   cualquier canal de 2 partículas — esto es consistente con cómo el propio paper enuncia las
   Conjeturas 1-3 ("para todo n≥3"), pero no se ha visto una declaración explícita de que el
   EXPONENTE del bound (`s²` específicamente, no otro) se mantenga igual para n≥5 — es la lectura
   más natural, pero vale la pena señalar el supuesto.
3. Dada la limpieza y robustez del resultado numérico (4 configuraciones independientes, mismo
   exponente a 3+ cifras), la conclusión de que el candidato de grado 9 viola el bound (sea `s²` u
   otro exponente par razonable) es sólida — lo que podría discutirse es únicamente el valor preciso
   del bound, no si el candidato lo excede.

## Archivos

- `scripts/regge_limit_n5.py` — construcción completa: clasificación de canales, referencia de
  gauge, evaluación del candidato de grado 9 y ajuste de pendiente.
