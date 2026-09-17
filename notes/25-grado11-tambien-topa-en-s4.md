# Fase 2 — Paso 24: buscando cancelación en grado 11 — corrección conceptual sobre cómo combinar
# grados, y hallazgo de que el piso de crecimiento sigue siendo s⁴

## Corrección conceptual antes de empezar (importante)

La idea original ("combinar el invariante de grado 9 con los de grado 11 para cancelar el exceso de
crecimiento") **está mal planteada**: 2 invariantes de grado TOTAL distinto en momentos no pueden
cancelarse entre sí en el límite de Regge — crecen como potencias DISTINTAS de `Λ` (`Λ⁹` vs `Λ¹¹`),
y la de mayor grado domina siempre a suficiente `Λ`, sin importar el coeficiente relativo (salvo que
ese coeficiente sea exactamente cero). La cancelación de un término dominante solo puede ocurrir
**entre invariantes del MISMO grado total**, cuando el espacio en ese grado tiene dimensión >1 —
exactamente el mecanismo que usa el paper de n=4 (sección 2.10): a 6 derivadas, la combinación
especial `stu·|e_S⟩` crece más lento que otros descendientes genéricos del mismo grado, porque son
combinaciones DENTRO del mismo nivel, no entre niveles distintos.

**Consecuencia inmediata para grado 9**: como su espacio de invariantes tiene dimensión 1 (parte 20),
no hay ningún otro elemento del MISMO grado con el que cancelar su crecimiento `s⁴` — su único
invariante **queda excluido sin remedio** de cualquier S-matrix física (su coeficiente en la
Lagrangiana/S-matrix debe ser exactamente cero).

**La pregunta que sí tiene sentido**: grado 11 tiene dimensión 3 (parte 21) — ¿alguna combinación
lineal DENTRO de ese espacio de 3 dimensiones crece más lento que `s⁴`?

## Método

Se identificaron 3 semillas explícitas que forman una base del espacio de grado 11 (extraídas como
las columnas pivote de la eliminación gaussiana de alta precisión ya usada en la parte 21):

```
B1 = F^i:F^j · F^l:F^m · C^k_{ij} · d(i,j)·d(i,l)
B2 = F^i:F^j · F^l:F^m · C^k_{ij} · d(i,l)²
B3 = F^i:F^j · F^l:F^m · C^k_{ij} · d(i,l)·d(i,m)
```

(todas simetrizadas sobre la órbita completa de S₅, misma técnica de siempre), evaluadas sobre el
límite de Regge de n=5 ya construido y validado (parte 23).

## Resultado

- **B1 y B2 crecen ambos como `E^10 ~ s^5`** (verificado en 3 configuraciones de parámetros
  distintas, slopes 9.999-10.001).
- **B3 crece como `E^8 ~ s^4`** por sí solo — de entrada, sin necesitar ninguna combinación, ya
  iguala el piso encontrado en grado 9.
- Se buscó la combinación `B1 + r·B2` que cancela el término líder `E^10` (razón `r` obtenida por
  extrapolación de Richardson sobre el cociente `B1/B2` a `E` grande, convergiendo establemente a
  `r≈0.6649`). **Resultado: la cancelación funciona** (el término `E^10` se suprime en 5 órdenes de
  magnitud), pero lo que queda **crece como `E^9 ~ s^4.5`** — sigue siendo PEOR que `B3` solo
  (`E^8~s^4`).

**Argumento de por qué `s⁴` es el piso genuino de este espacio de 3 dimensiones** (no solo lo que se
encontró probando, sino por qué no puede haber nada mejor): cualquier combinación
`c₁B1+c₂B2+c₃B3` que no sea pura `B3` (es decir, con `c₁` o `c₂` no nulos) necesariamente contiene
el término `E^10` de `B1,B2`, salvo que se ajuste `c₂/c₁` exactamente a la razón de cancelación — pero
esa combinación deja un residuo `E^9` que NO depende de `c₃` (viene enteramente de la parte subléder
de `B1,B2`, que `B3` — de grado E^8 — no puede tocar). Por lo tanto, **cualquier combinación con
`c₁≠0` crece como mínimo `E^9`, peor que `E^8`**. La única forma de alcanzar el piso `E^8` es poner
`c₁=c₂=0`, es decir, **usar B3 solo**. `s⁴` es entonces el mínimo genuino alcanzable en todo el
espacio de grado 11, no una limitación del método de búsqueda.

## Conclusión: grados 9 y 11 topan igual, en `s⁴` — ninguno satisface Regge growth

**Tanto el único invariante de grado 9 como el mejor invariante alcanzable en grado 11 (dimensión 3)
crecen exactamente como `s⁴`**, excediendo el bound `s²` por el mismo margen. Esto es un hallazgo
más fuerte que "grado 9 falla": sugiere que **`s⁴` podría ser un piso estructural de esta familia de
construcciones** (basadas en los átomos `F^i:F^j` y `C^k_{ab}`), no un accidente del grado más bajo.

**Hipótesis a investigar** (no confirmada, especulativa): la presencia del building block
`C^k_{ab}` (el "singlete", lineal en una sola polarización) podría ser estructuralmente responsable
de este piso — es el único átomo que no es manifiestamente "grado 2 por pareja", y quizás cualquier
construcción que lo use hereda un crecimiento mínimo de `s⁴` en este canal. Si esto se confirmara,
implicaría que ninguna S-matrix construida con la topología "parejas + singlete" puede satisfacer
Regge growth para n=5, y haría falta una estructura fundamentalmente distinta (quizás el objeto
"sandwich" `p_a.F^i.F^j.p_b` mencionado como pendiente en partes anteriores, u otra cosa).

## Qué falta (próximo paso concreto)

1. Verificar la hipótesis: evaluar el crecimiento del objeto "sandwich" (`p_a.F^i.F^j.p_b`, nunca
   explorado numéricamente en esta sesión) en este mismo límite de Regge, para ver si evita el piso
   de `s⁴` o lo hereda también.
2. Si la hipótesis se confirma (todo lo construido con `C^k_{ab}` topa en `s⁴`), la conclusión
   parcial honesta sería: **con la topología "2 parejas + 1 singlete" (la única explorada hasta
   ahora, y la de grado mínimo conocido, parte 19), no existe ningún candidato de fotones n=5 que
   satisfaga Regge growth** — habría que buscar topologías genuinamente distintas.
3. Si aparece un candidato de grado más alto (13, 15, ...) con dimensión >1 que sí logre bajar de
   `s⁴`, seguir esa línea en cambio.

## Advertencia honesta

Este resultado depende de: (a) la validez de la generalización del límite de Regge a n=5 ya
construida (parte 23, validada contra n=4 pero es una elección propia), (b) que el bound relevante
sea efectivamente `s²` para n=5 (mismo exponente que n=4, supuesto natural pero no confirmado en la
literatura). El **argumento de que `s⁴` es el piso DENTRO del espacio de grado 11** (dado el límite
de Regge y las 3 semillas encontradas) es sólido y no depende de esos supuestos externos — lo que
podría cambiar con otro límite de Regge es el VALOR del exponente medido, no la conclusión relativa
de que grado 9 y grado 11 topan en el mismo lugar.

## Archivos

- Cálculos hechos inline (reutilizando `regge_limit_n5.py` y las semillas identificadas via
  `hilbert_series_grado11_altaprecision.py`); no se generó un script nuevo separado dado que son
  extensiones cortas de scripts ya existentes.
