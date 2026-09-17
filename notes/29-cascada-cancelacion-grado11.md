# Fase 2 — Paso 29: grado 11 muestra la misma "cascada" de cancelación que grado 9 — bajado de
# s⁵ a ~s^3.5, con el último paso pendiente de precisión

## Resultado

Se evaluaron los 15 elementos base de grado 11 (parte 28) en el límite de Regge de n=5:
**13 de los 15 crecen como `s⁵`, y 2 (S3, S7) ya de entrada crecen como `s⁴`** (el nivel encontrado
en grado 9).

**Se encontraron relaciones exactas entre coeficientes líderes** (análogas al `A1=2A2` de grado 9):
`A(S8) = -0.5·A(S1)` y `A(S14) = A(S15)` (coincidencia limpia a 6+ cifras). Esto da 2 combinaciones
que cancelan el término `s⁵` de inmediato:

```
w1 := S1 + 2·S8   ~ s⁴  (pendiente 3.9999, limpio)
w2 := S14 - S15   ~ s⁴  (pendiente 4.0001, limpio)
```

**Segundo paso**: con 4 elementos disponibles en el nivel `s⁴` (S3, S7, w1, w2), se buscó
cancelación adicional. La combinación `w1 + c·S3` (con `c` obtenido por extrapolación de Richardson
de la razón `w1/S3`, convergiendo a `c≈-30.905`) **baja el crecimiento a `~s^3.48`** — mejora real
(de `s⁴` a algo cercano a `s^3.5`, el mismo nivel intermedio encontrado en grado 9), pero **la
pendiente no es tan limpia como en los pasos anteriores** (6.965 en vez de un entero exacto),
señal de que el coeficiente `c` todavía no está calculado con precisión suficiente (mismo problema
ya documentado en la parte 27 para la 2ª dirección de grado 9).

## Patrón que emerge (para ambos grados 9 y 11)

Existe una estructura de "cascada" clara: `s⁵ → s⁴ → ~s^3.5 → ?`, donde cada paso de cancelación es
más difícil de calcular con precisión que el anterior (los coeficientes de cancelación se vuelven
más sensibles/irracionales, no combinaciones simples como `-2` o `0.5`). Esto es consistente con
cómo funciona la teoría de invariantes en general: mientras más "profunda" la cancelación, más fina
la relación algebraica necesaria — y confirma que **hace falta una derivación simbólica exacta**
(no más ajustes numéricos con Richardson) para llegar con confianza hasta el final (¿existe un
candidato que llegue a `s²`, el bound real?).

## Estado honesto

- **Confirmado con rigor numérico sólido**: grado 11 baja de `s⁵` a `s⁴` (2 combinaciones limpias,
  pendientes exactas).
- **Sugerido, no confirmado con rigor suficiente**: una bajada adicional hasta `~s^3.5` — la
  pendiente medida (6.97 en vez de 7.00) indica que el coeficiente usado no es exacto.
- **No se ha llegado a `s²`** en ninguno de los 2 grados explorados (9 y 11) — sigue siendo la meta
  pendiente.

## Recomendación para la continuidad del proyecto

Dado que la vía puramente numérica está llegando a un punto de rendimientos decrecientes (cada paso
de cancelación requiere más precisión que la anterior, y el ajuste por Richardson numérico no
escala bien más allá de 2-3 niveles), **el paso más valioso ahora es intentar una derivación
simbólica** de los coeficientes de cancelación — usando el anillo de Mandelstam ya calculado (partes
3-8) y las fórmulas cerradas ya derivadas para F^i:F^j (partes 16-17) y el sandwich (parte 25) de
forma exacta, en vez de numérica. Esto es un cambio de metodología significativo, probablemente
mejor abordado como el foco de una sesión dedicada, no como una extensión más del enfoque numérico
actual.

## Archivos

- `scripts/grado11_regge_cancelacion.py` — los 15 elementos base evaluados en el límite de Regge,
  con las combinaciones `w1, w2` y el intento de segunda cancelación.
