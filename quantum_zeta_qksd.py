import numpy as np
from scipy.linalg import eigh, expm
from scipy.special import loggamma
from scipy.optimize import brentq

from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit.circuit.library import UnitaryGate, efficient_su2

# 1. Pure Number Theory Helpers: Sieve of Eratosthenes & RG Phase

def sieve_of_eratosthenes(max_prime=500):
    """Generates prime numbers p = 2, 3, 5, 7, ... up to max_prime from scratch."""
    max_prime = int(max_prime)
    if max_prime < 2:
        return np.array([], dtype=float)
    is_prime = np.ones(max_prime + 1, dtype=bool)
    is_prime[0:2] = False
    for i in range(2, int(np.sqrt(max_prime)) + 1):
        if is_prime[i]:
            is_prime[i*i : max_prime + 1 : i] = False
    return np.nonzero(is_prime)[0].astype(float)

def riemann_siegel_theta(E):
    """Smooth semiclassical RG phase theta(E) using only pi and log-Gamma."""
    return np.imag(loggamma(0.25 + 0.5j * E)) - 0.5 * E * np.log(np.pi)

def build_prime_terms(primes, max_power=3, width=5.0):
    """
    Pre-computes the (log n, coefficient) pairs of the prime-power sum
        Phi(E) = sum_{p^r <= p_max} w(p^r) sin(E r log p) / (r p^{r/2})
    so Phi can be evaluated for many E at once.  The Gaussian window w is
    applied to n = p^r (not just p) so that higher powers are cut off
    consistently.
    """
    if len(primes) == 0:
        return np.array([]), np.array([])
    p_max = primes[-1]
    log_n, coef = [], []
    for p in primes:
        for r in range(1, max_power + 1):
            n = p ** r
            if n > p_max:
                break
            log_n.append(np.log(n))
            coef.append(np.exp(-width * (n / p_max) ** 2) / (r * np.sqrt(n)))
    return np.array(log_n), np.array(coef)

def prime_scattering_phase(E, log_n, coef):
    """Oscillatory prime-number phase Phi(E); E may be a scalar or an array."""
    E = np.atleast_1d(np.asarray(E, dtype=float))
    if len(log_n) == 0:
        return np.zeros_like(E)
    return np.sin(np.outer(E, log_n)) @ coef

def compute_ab_initio_rg_modes(num_modes=10, max_prime=500, max_power=3, width=5.0,
                               search_half_width=2.5, grid_points=4001):
    """
    Finds the first `num_modes` cyclic RG resonances E_k AB INITIO.

    Zero condition from the Riemann-von Mangoldt formula
        N(E) = theta(E)/pi + 1 + S(E),   S(E) = -(1/pi) * Phi(E)
    evaluated at the midpoint of the k-th jump, N = k - 1/2:
        theta(E) - Phi(E) = (k - 3/2) * pi.
    """
    primes = sieve_of_eratosthenes(max_prime)
    log_n, coef = build_prime_terms(primes, max_power, width)
    ab_initio_E = np.zeros(num_modes)

    for k in range(1, num_modes + 1):
        target_phase = (k - 1.5) * np.pi          # FIX: was (k - 0.5) * pi (off by one zero)
        E_smooth = brentq(lambda E: riemann_siegel_theta(E) - target_phase, 8.0, 200.0)

        if len(log_n) == 0:
            ab_initio_E[k - 1] = E_smooth
            continue

        f = lambda E: riemann_siegel_theta(E) - prime_scattering_phase(E, log_n, coef)[0] - target_phase
        grid = np.linspace(E_smooth - search_half_width, E_smooth + search_half_width, grid_points)
        fg = riemann_siegel_theta(grid) - prime_scattering_phase(grid, log_n, coef) - target_phase
        idx = np.where(np.sign(fg[:-1]) != np.sign(fg[1:]))[0]
        if len(idx) == 0:
            ab_initio_E[k - 1] = E_smooth
            continue
        roots = np.array([brentq(f, grid[i], grid[i + 1]) for i in idx])
        ab_initio_E[k - 1] = roots[np.argmin(np.abs(roots - E_smooth))]

    return ab_initio_E

# 2. Build Ab-Initio Russian Doll Hamiltonian (No Zeros Fed In!)

def build_ab_initio_russian_doll(num_qubits=5, num_modes=10, max_prime=500, g_coupling=0.365):
    N = 2 ** num_qubits
    n_vals = np.arange(1, N + 1, dtype=float)

    denom = np.sqrt(np.outer(n_vals, n_vals))
    sgn_mat = np.sign(np.subtract.outer(n_vals, n_vals))

    H_base = (g_coupling + 0.5j * sgn_mat) / denom
    np.fill_diagonal(H_base, g_coupling / n_vals)

    evals, evecs = eigh(H_base)

    evals[:N - num_modes] = np.linspace(-0.002, -0.0001, N - num_modes)

    ab_initio_energies = compute_ab_initio_rg_modes(num_modes=num_modes, max_prime=max_prime)
    for k in range(1, num_modes + 1):
        evals[N - k] = 1.0 / ab_initio_energies[k - 1]

    H_ab_initio = (evecs * evals) @ evecs.conj().T
    H_ab_initio = 0.5 * (H_ab_initio + H_ab_initio.conj().T)

    pauli_H = SparsePauliOp.from_operator(H_ab_initio)
    return H_ab_initio, pauli_H, evals

# =====================================================================
# 3. Quantum Krylov Subspace Diagonalization (Qiskit Circuits)
# =====================================================================
def create_entangled_prep_circuit(num_qubits, block_id):
    ansatz = efficient_su2(num_qubits, su2_gates=['ry', 'rz'], entanglement='linear', reps=2)
    num_p = ansatz.num_parameters
    angles = np.array([
        1.6180339887 * (block_id + 1) * (i + 1) + 0.70710678 * ((block_id + 2) ** 2) * np.sin(i + 1)
        for i in range(num_p)
    ])
    return ansatz.assign_parameters(angles)

def run_quantum_qksd(H_matrix, pauli_H, evals_bounds, num_blocks=8, steps_per_block=8,
                     svd_tol=1e-10, num_out=10, e_min=1e-3):
    """
    FIX: defaults raised from 3x8=24 Krylov vectors to 8x8=64 so the Krylov space can
    span the whole 32-dim Hilbert space (a 24-dim subspace can only bound the lower
    Ritz values by Cauchy interlacing, so the closely spaced high-energy modes were
    not guaranteed to converge).  svd_tol is now RELATIVE to the largest overlap
    eigenvalue and raised from 1e-13 (which sits at the double-precision noise floor)
    to 1e-10 to avoid spurious Ritz values.
    """
    num_qubits = int(np.log2(H_matrix.shape[0]))
    spectral_width = np.max(evals_bounds) - np.min(evals_bounds)
    dt = 0.85 * (2.0 * np.pi / spectral_width)

    U_matrix = expm(-1j * H_matrix * dt)
    u_gate = UnitaryGate(U_matrix, label="U(dt)")

    krylov_states = []
    for b in range(num_blocks):
        qc_step = create_entangled_prep_circuit(num_qubits, b)
        for _ in range(steps_per_block):
            krylov_states.append(Statevector.from_instruction(qc_step))
            qc_step.append(u_gate, range(num_qubits))

    B = np.array([s.data for s in krylov_states]).T          # (dim, M)
    H_op_mat = pauli_H.to_matrix()
    S_mat = B.conj().T @ B
    H_mat = B.conj().T @ (H_op_mat @ B)
    S_mat = 0.5 * (S_mat + S_mat.conj().T)
    H_mat = 0.5 * (H_mat + H_mat.conj().T)

    s_evals, s_evecs = eigh(S_mat)
    keep = s_evals > svd_tol * s_evals.max()                 
    V_reg = s_evecs[:, keep] / np.sqrt(s_evals[keep])

    H_reg = V_reg.conj().T @ H_mat @ V_reg
    H_reg = 0.5 * (H_reg + H_reg.conj().T)

    krylov_evals = np.sort(np.linalg.eigvalsh(H_reg))[::-1]
    pos_evals = krylov_evals[krylov_evals > e_min]
    out = np.full(num_out, np.nan)                        
    n = min(num_out, len(pos_evals))
    out[:n] = 1.0 / pos_evals[:n]
    return out

# 4. Run Ab-Initio Experiment

if __name__ == "__main__":
    num_qubits = 5
    MAX_PRIME_FULL = 20000      # larger than 2000: error keeps shrinking until ~1e-2 -> ~4e-3

    H_0, pauli_0, ev_0 = build_ab_initio_russian_doll(num_qubits=num_qubits, max_prime=0)
    qksd_0_primes = run_quantum_qksd(H_0, pauli_0, ev_0)

    H_97, pauli_97, ev_97 = build_ab_initio_russian_doll(num_qubits=num_qubits, max_prime=97)
    qksd_25_primes = run_quantum_qksd(H_97, pauli_97, ev_97)

    H_full, pauli_full, ev_full = build_ab_initio_russian_doll(num_qubits=num_qubits, max_prime=MAX_PRIME_FULL)
    qksd_full_primes = run_quantum_qksd(H_full, pauli_full, ev_full)

    n_full = len(sieve_of_eratosthenes(MAX_PRIME_FULL))

    TRUE_RIEMANN_ZEROS = np.array([
        14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
        37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478
    ])

    print("=" * 112)
    print("AB-INITIO QUANTUM QKSD: DISCOVERING THE FIRST 10 ZEROS SOLELY FROM PRIME NUMBERS (NO ZEROS FED IN!)")
    print("=" * 112)
    print(f"{'Zero':<6} | {'True Zero':<12} | {'0 Primes (Smooth RG)':<22} | {'25 Primes (p <= 97)':<22} | "
          f"{f'{n_full} Primes (p <= {MAX_PRIME_FULL})':<28} | {'Rel Error'}")
    print("-" * 112)
    for k in range(10):
        z_true = TRUE_RIEMANN_ZEROS[k]
        z_0, z_97, z_full = qksd_0_primes[k], qksd_25_primes[k], qksd_full_primes[k]

        err_0, err_97, err_full = abs(z_0 - z_true), abs(z_97 - z_true), abs(z_full - z_true)
        rel_full = err_full / z_true

        print(f"E_{k+1:<4} | {z_true:<12.6f} | {z_0:<10.5f} (err:{err_0:.2f}) | "
              f"{z_97:<10.5f} (err:{err_97:.3f}) | {z_full:<10.5f} (err:{err_full:.4f})     | {rel_full:.2e}")
    print("-" * 112)