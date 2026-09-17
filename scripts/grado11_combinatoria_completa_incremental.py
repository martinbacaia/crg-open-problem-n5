"""
Version incremental (resumible) de grado11_combinatoria_completa.py -- se descubrio que dejar esto
corriendo con run_in_background (proceso totalmente desatendido) es poco confiable en este entorno
(se atasca sin usar CPU durante largos ratos, por razones de scheduling del sistema, no del codigo
-- verificado: la MISMA computacion corrida en foreground con timeout explicito SI funciona limpio,
sin colgarse, con tiempos por muestra de 50-135s). Por eso esta version guarda resultados
incrementales en un pickle y se corre en VARIAS llamadas foreground cortas (cada una agrega unas
pocas muestras), nunca en background desatendido.

Uso: python grado11_combinatoria_completa_incremental.py <n_muestras_a_agregar>
Guarda/actualiza "combinatoria_completa_rows.pkl" con las filas ya calculadas.
"""
import sys, time, pickle, os
from grado11_combinatoria_completa import build_seeds
from grado11_con_familia_D import gen_kinematics, orbit_sum, mp_matrix_rank

PKL = "combinatoria_completa_rows.pkl"

def load():
    if os.path.exists(PKL):
        with open(PKL, "rb") as f:
            return pickle.load(f)
    return {"seed_names": None, "rows": {}}  # rows: {kseed: row}

def save(data):
    with open(PKL, "wb") as f:
        pickle.dump(data, f)

if __name__ == "__main__":
    n_to_add = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    seeds = build_seeds()
    seed_names = list(seeds.keys())

    data = load()
    if data["seed_names"] is None:
        data["seed_names"] = seed_names
    assert data["seed_names"] == seed_names, "seed set changed, no se puede resumir"

    done_seeds = sorted(data["rows"].keys())
    next_kseed = (done_seeds[-1] + 1) if done_seeds else 1
    print(f"ya calculadas: {len(done_seeds)} muestras. Agregando {n_to_add} mas, desde kseed={next_kseed}", flush=True)

    for kseed in range(next_kseed, next_kseed + n_to_add):
        t0 = time.time()
        ps, eps, eta = gen_kinematics(kseed)
        row = [orbit_sum(eta, ps, eps, seeds[name]) for name in seed_names]
        data["rows"][kseed] = row
        save(data)
        print(f"  muestra {kseed} calculada y guardada ({time.time()-t0:.1f}s)", flush=True)

    n_total = len(data["rows"])
    print(f"total acumulado: {n_total} muestras", flush=True)

    if len(sys.argv) > 2 and sys.argv[2] == "rank":
        rows = [data["rows"][k] for k in sorted(data["rows"].keys())]
        norms = [mp_matrix_rank.__globals__["mp"].sqrt(sum(rows[r][c] ** 2 for r in range(n_total))) for c in range(len(seed_names))]
        mp = mp_matrix_rank.__globals__["mp"]
        nz_idx = [c for c, n in enumerate(norms) if n > mp.mpf('1e-25')]
        print(f"Semillas no nulas: {len(nz_idx)} de {len(seed_names)}")
        Mn = [[rows[r][c] / norms[c] for c in nz_idx] for r in range(n_total)]
        rank = mp_matrix_rank(Mn, verbose=True)
        print(f"\nRango (combinatoria completa, {n_total} muestras) = {rank}")
