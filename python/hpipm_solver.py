#!/usr/bin/env python3
"""
HPIPM / cvxpy QCQP Solver Interface

Solves QCQP problems arising from RSS controller.
Supports two backends:
  1. HPIPM (hpipm_python with compiled C library) - high performance
  2. cvxpy (with ECOS/Clarabel/OSQP) - portable fallback

QCQP problem form:
  min  0.5*x'*H*x + g'*x
  s.t. A*x = b                  (linear equality)
       C*x <= d                 (linear inequality)
       lb <= x <= ub            (box bounds, optional)
       0.5*x'*Hq_i*x + gq_i'*x <= uq_i  (quadratic inequality, optional)

Called from MATLAB via system():
  python hpipm_solver.py --input qp_input.json --output qp_output.json
"""

import numpy as np
import sys
import json
import argparse
import time
from pathlib import Path

# --- Solver availability ---
HPIPM_AVAILABLE = False
CVXPY_AVAILABLE = False

try:
    from hpipm_python import *
    from hpipm_python.common import *
    HPIPM_AVAILABLE = True
except ImportError:
    pass

try:
    import cvxpy as cp
    CVXPY_AVAILABLE = True
except ImportError:
    pass


def solve_qcqp(H, g, A=None, b=None, C=None, d=None,
               lb=None, ub=None, Hq=None, gq=None, uq=None,
               solver='auto', verbose=False):
    """
    Solve a QCQP problem.

    Parameters
    ----------
    H : ndarray (n, n)  - Hessian (must be symmetric PSD)
    g : ndarray (n,)    - linear cost term
    A : ndarray (me, n) or None - equality constraint matrix
    b : ndarray (me,)   or None - equality constraint RHS
    C : ndarray (mg, n) or None - inequality constraint matrix (C*x <= d)
    d : ndarray (mg,)   or None - inequality constraint RHS
    lb : ndarray (n,)   or None - lower bounds
    ub : ndarray (n,)   or None - upper bounds
    Hq : list[ndarray]  or None - Hessian matrices for quad constraints
    gq : list[ndarray]  or None - linear terms for quad constraints
    uq : list[float]    or None - upper bounds for quad constraints
    solver : str - 'auto', 'hpipm', or 'cvxpy'
    verbose : bool

    Returns
    -------
    dict with keys: x, status, obj_value, solve_time, solver_used
    """
    H = np.asarray(H, dtype=np.float64)
    g = np.asarray(g, dtype=np.float64).flatten()
    n = H.shape[0]

    # Ensure H is symmetric
    H = 0.5 * (H + H.T)

    # Normalize empty constraints
    nq = 0
    if Hq is not None and len(Hq) > 0:
        nq = len(Hq)
    has_eq = A is not None and b is not None and len(b) > 0
    has_ineq = C is not None and d is not None and len(d) > 0

    # Select solver — track requested vs actual for audit
    requested_solver = solver
    fallback_used = False

    if solver == 'auto':
        if nq > 0:
            solver = 'hpipm' if HPIPM_AVAILABLE else 'cvxpy'
        else:
            solver = 'hpipm' if HPIPM_AVAILABLE else 'cvxpy'

    if solver == 'hpipm' and not HPIPM_AVAILABLE:
        print("HPIPM not available, falling back to cvxpy", file=sys.stderr)
        solver = 'cvxpy'
        fallback_used = True

    if solver == 'cvxpy' and not CVXPY_AVAILABLE:
        raise RuntimeError("Neither HPIPM nor cvxpy is available. Install cvxpy: pip install cvxpy")

    # Dispatch
    if solver == 'hpipm':
        if nq > 0:
            result = _solve_qcqp_hpipm(H, g, A, b, C, d, lb, ub, Hq, gq, uq, verbose)
        else:
            result = _solve_qp_hpipm(H, g, A, b, C, d, lb, ub, verbose)
    else:
        result = _solve_qcqp_cvxpy(H, g, A, b, C, d, lb, ub, Hq, gq, uq, verbose)

    # Add audit metadata to every result
    result['requestedSolver'] = requested_solver
    result['actualSolver'] = solver
    result['fallbackUsed'] = fallback_used
    return result


# =====================================================================
# HPIPM dense QP solver
# =====================================================================

def _solve_qp_hpipm(H, g, A, b, C, d, lb, ub, verbose=False):
    """Solve QP using HPIPM dense_qp interface."""
    n = H.shape[0]
    ne = A.shape[0] if A is not None else 0
    ng = C.shape[0] if C is not None else 0
    nb = 0
    if lb is not None or ub is not None:
        nb = n  # all variables have bounds

    # Dimension
    dim = hpipm_dense_qp_dim()
    dim.set('nv', n)
    dim.set('ne', ne)
    dim.set('nb', nb)
    dim.set('ng', ng)

    # QP data
    qp = hpipm_dense_qp(dim)
    qp.set('H', H)
    qp.set('g', g)

    if ne > 0:
        qp.set('A', np.asarray(A, dtype=np.float64))
        qp.set('b', np.asarray(b, dtype=np.float64))

    if ng > 0:
        C_arr = np.asarray(C, dtype=np.float64)
        d_arr = np.asarray(d, dtype=np.float64)
        qp.set('C', C_arr)
        qp.set('ug', d_arr)
        qp.set('lg', -np.inf * np.ones(ng))

    if nb > 0:
        idxb = np.arange(n)
        qp.set('idxb', idxb)
        if lb is not None:
            qp.set('lb', np.asarray(lb, dtype=np.float64))
        else:
            qp.set('lb', -np.inf * np.ones(n))
        if ub is not None:
            qp.set('ub', np.asarray(ub, dtype=np.float64))
        else:
            qp.set('ub', np.inf * np.ones(n))

    # Solution
    qp_sol = hpipm_dense_qp_sol(dim)

    # Solver arguments — use 'balance' mode for better robustness
    arg = hpipm_dense_qp_solver_arg(dim, 'balance')
    arg.set('iter_max', 300)
    arg.set('tol_stat', 1e-6)
    arg.set('tol_eq', 1e-6)
    arg.set('tol_ineq', 1e-6)
    arg.set('tol_comp', 1e-6)

    # Solve
    solver = hpipm_dense_qp_solver(dim, arg)
    t0 = time.time()
    solver.solve(qp, qp_sol)
    solve_time = time.time() - t0

    # Extract results
    x = qp_sol.get('v').flatten()
    status = solver.get('status')
    obj_value = 0.5 * x @ H @ x + g @ x

    status_str = 'Solved' if status == 0 else f'Failed({status})'

    # If balanced mode failed, retry with robust mode
    if status != 0:
        print(f"HPIPM QP balanced mode failed (status={status}), retrying with robust mode...", file=sys.stderr)
        arg2 = hpipm_dense_qp_solver_arg(dim, 'robust')
        arg2.set('iter_max', 500)
        arg2.set('tol_stat', 1e-6)
        arg2.set('tol_eq', 1e-6)
        arg2.set('tol_ineq', 1e-6)
        arg2.set('tol_comp', 1e-6)
        solver2 = hpipm_dense_qp_solver(dim, arg2)
        t1 = time.time()
        solver2.solve(qp, qp_sol)
        solve_time += time.time() - t1
        x = qp_sol.get('v').flatten()
        status = solver2.get('status')
        obj_value = 0.5 * x @ H @ x + g @ x
        status_str = 'Solved' if status == 0 else f'Failed({status})'

    return {
        'x': x,
        'status': status_str,
        'obj_value': float(obj_value),
        'solve_time': float(solve_time),
        'solver_used': 'hpipm_dense_qp'
    }


# =====================================================================
# HPIPM dense QCQP solver
# =====================================================================

def _solve_qcqp_hpipm(H, g, A, b, C, d, lb, ub, Hq, gq, uq, verbose=False):
    """Solve QCQP using HPIPM dense_qcqp interface."""
    n = H.shape[0]
    ne = A.shape[0] if A is not None else 0
    ng = C.shape[0] if C is not None else 0
    nq = len(Hq) if Hq is not None else 0
    nb = 0
    if lb is not None or ub is not None:
        nb = n

    # Dimension
    dim = hpipm_dense_qcqp_dim()
    dim.set('nv', n)
    dim.set('ne', ne)
    dim.set('nb', nb)
    dim.set('ng', ng)
    dim.set('nq', nq)

    # QCQP data
    qcqp = hpipm_dense_qcqp(dim)
    qcqp.set('H', H)
    qcqp.set('g', g)

    if ne > 0:
        qcqp.set('A', np.asarray(A, dtype=np.float64))
        qcqp.set('b', np.asarray(b, dtype=np.float64))

    if ng > 0:
        C_arr = np.asarray(C, dtype=np.float64)
        d_arr = np.asarray(d, dtype=np.float64)
        qcqp.set('C', C_arr)
        qcqp.set('ug', d_arr)
        qcqp.set('lg', -np.inf * np.ones(ng))

    if nb > 0:
        idxb = np.arange(n)
        qcqp.set('idxb', idxb)
        if lb is not None:
            qcqp.set('lb', np.asarray(lb, dtype=np.float64))
        else:
            qcqp.set('lb', -np.inf * np.ones(n))
        if ub is not None:
            qcqp.set('ub', np.asarray(ub, dtype=np.float64))
        else:
            qcqp.set('ub', np.inf * np.ones(n))

    # Quadratic constraints: stack Hq and gq
    # [Bug fix] gq_stacked must be (n, nq) matrix for dense QCQP interface,
    # not a 1D vector from hstack. Use column_stack instead.
    if nq > 0:
        Hq_stacked = np.hstack([np.asarray(hq, dtype=np.float64) for hq in Hq])
        gq_stacked = np.column_stack([
            np.asarray(gqi, dtype=np.float64).reshape(-1)
            for gqi in gq
        ])
        uq_arr = np.asarray(uq, dtype=np.float64).flatten()
        qcqp.set('Hq', Hq_stacked)
        qcqp.set('gq', gq_stacked)
        qcqp.set('uq', uq_arr)

    # Solution
    qcqp_sol = hpipm_dense_qcqp_sol(dim)

    # Solver arguments — use 'balance' mode for better robustness
    arg = hpipm_dense_qcqp_solver_arg(dim, 'balance')
    arg.set('iter_max', 300)
    arg.set('tol_stat', 1e-6)
    arg.set('tol_eq', 1e-6)
    arg.set('tol_ineq', 1e-6)
    arg.set('tol_comp', 1e-6)

    # Solve
    solver = hpipm_dense_qcqp_solver(dim, arg)
    t0 = time.time()
    solver.solve(qcqp, qcqp_sol)
    solve_time = time.time() - t0

    # Extract results
    x = qcqp_sol.get('v').flatten()
    status = solver.get('status')
    obj_value = 0.5 * x @ H @ x + g @ x

    status_str = 'Solved' if status == 0 else f'Failed({status})'

    # If balanced mode failed, retry with robust mode
    if status != 0:
        print(f"HPIPM QCQP balanced mode failed (status={status}), retrying with robust mode...", file=sys.stderr)
        arg2 = hpipm_dense_qcqp_solver_arg(dim, 'robust')
        arg2.set('iter_max', 500)
        arg2.set('tol_stat', 1e-6)
        arg2.set('tol_eq', 1e-6)
        arg2.set('tol_ineq', 1e-6)
        arg2.set('tol_comp', 1e-6)
        solver2 = hpipm_dense_qcqp_solver(dim, arg2)
        t1 = time.time()
        solver2.solve(qcqp, qcqp_sol)
        solve_time += time.time() - t1
        x = qcqp_sol.get('v').flatten()
        status = solver2.get('status')
        obj_value = 0.5 * x @ H @ x + g @ x
        status_str = 'Solved' if status == 0 else f'Failed({status})'

    return {
        'x': x,
        'status': status_str,
        'obj_value': float(obj_value),
        'solve_time': float(solve_time),
        'solver_used': 'hpipm_dense_qcqp'
    }


# =====================================================================
# cvxpy solver (fallback)
# =====================================================================

def _solve_qcqp_cvxpy(H, g, A, b, C, d, lb, ub, Hq, gq, uq, verbose=False):
    """Solve QCQP using cvxpy."""
    n = H.shape[0]

    x = cp.Variable(n)

    # Objective
    objective = cp.Minimize(0.5 * cp.quad_form(x, H, assume_PSD=True) + g @ x)

    constraints = []

    # Equality constraints
    if A is not None and b is not None and len(b) > 0:
        constraints.append(np.asarray(A) @ x == np.asarray(b).flatten())

    # Linear inequality constraints
    if C is not None and d is not None and len(d) > 0:
        constraints.append(np.asarray(C) @ x <= np.asarray(d).flatten())

    # Box constraints
    if lb is not None:
        constraints.append(x >= np.asarray(lb).flatten())
    if ub is not None:
        constraints.append(x <= np.asarray(ub).flatten())

    # Quadratic inequality constraints
    nq = 0
    if Hq is not None and len(Hq) > 0:
        nq = len(Hq)
        for i in range(nq):
            Hq_i = np.asarray(Hq[i], dtype=np.float64)
            Hq_i = 0.5 * (Hq_i + Hq_i.T)  # ensure symmetric
            gq_i = np.asarray(gq[i], dtype=np.float64).flatten()
            constraints.append(
                0.5 * cp.quad_form(x, Hq_i, assume_PSD=True) + gq_i @ x <= float(uq[i])
            )

    # Solve
    prob = cp.Problem(objective, constraints)

    t0 = time.time()
    try:
        # Try Clarabel first (good for QCQP), then ECOS, then SCS
        prob.solve(solver=cp.CLARABEL, verbose=verbose)
    except cp.SolverError:
        try:
            prob.solve(solver=cp.ECOS, verbose=verbose)
        except cp.SolverError:
            prob.solve(solver=cp.SCS, verbose=verbose)
    solve_time = time.time() - t0

    if x.value is None:
        return {
            'x': np.full(n, np.nan),
            'status': 'Failed',
            'obj_value': float('inf'),
            'solve_time': float(solve_time),
            'solver_used': 'cvxpy'
        }

    x_val = x.value.flatten()
    obj_value = 0.5 * x_val @ H @ x_val + g @ x_val

    # Map cvxpy status
    status_map = {
        cp.OPTIMAL: 'Solved',
        cp.OPTIMAL_INACCURATE: 'Inaccurate/Solved',
        cp.INFEASIBLE: 'Infeasible',
        cp.INFEASIBLE_INACCURATE: 'Inaccurate/Infeasible',
        cp.UNBOUNDED: 'Unbounded',
    }
    status_str = status_map.get(prob.status, str(prob.status))

    return {
        'x': x_val,
        'status': status_str,
        'obj_value': float(obj_value),
        'solve_time': float(solve_time),
        'solver_used': 'cvxpy'
    }


# =====================================================================
# JSON I/O for MATLAB interface
# =====================================================================

def load_qcqp_from_json(json_file):
    """Load QCQP problem from JSON file (written by MATLAB)."""
    with open(json_file, 'r') as f:
        data = json.load(f)

    H = np.array(data['H'], dtype=np.float64)
    g = np.array(data['g'], dtype=np.float64).flatten()

    A = np.array(data.get('A', []), dtype=np.float64) if data.get('A') and len(data['A']) > 0 else None
    b = np.array(data.get('b', []), dtype=np.float64).flatten() if data.get('b') and len(data['b']) > 0 else None
    C = np.array(data.get('C', []), dtype=np.float64) if data.get('C') and len(data['C']) > 0 else None
    d = np.array(data.get('d', []), dtype=np.float64).flatten() if data.get('d') and len(data['d']) > 0 else None
    lb = np.array(data.get('lb', []), dtype=np.float64).flatten() if data.get('lb') and len(data.get('lb', [])) > 0 else None
    ub = np.array(data.get('ub', []), dtype=np.float64).flatten() if data.get('ub') and len(data.get('ub', [])) > 0 else None

    # Quadratic constraints
    Hq_list = data.get('Hq', [])
    gq_list = data.get('gq', [])
    uq_list = data.get('uq', [])

    Hq = [np.array(h, dtype=np.float64) for h in Hq_list] if Hq_list and len(Hq_list) > 0 else None
    gq = [np.array(g, dtype=np.float64).flatten() for g in gq_list] if gq_list and len(gq_list) > 0 else None
    uq = [float(u) for u in uq_list] if uq_list and len(uq_list) > 0 else None

    return H, g, A, b, C, d, lb, ub, Hq, gq, uq


def save_solution_to_json(solution, output_file):
    """Save solution to JSON file for MATLAB to read."""
    output_data = {
        'x': solution['x'].tolist(),
        'status': solution['status'],
        'obj_value': float(solution['obj_value']),
        'solve_time': float(solution['solve_time']),
        'solver_used': solution['solver_used'],
        'requestedSolver': solution.get('requestedSolver', 'auto'),
        'actualSolver': solution.get('actualSolver', solution['solver_used']),
        'fallbackUsed': solution.get('fallbackUsed', False),
    }

    with open(output_file, 'w') as f:
        json.dump(output_data, f)


# =====================================================================
# CLI entry point
# =====================================================================

def main():
    parser = argparse.ArgumentParser(description='QCQP Solver (HPIPM/cvxpy)')
    parser.add_argument('--input', required=True, help='Input JSON file')
    parser.add_argument('--output', required=True, help='Output JSON file')
    parser.add_argument('--solver', default='auto', choices=['auto', 'hpipm', 'cvxpy'])
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()

    # Load problem
    try:
        H, g, A, b, C, d, lb, ub, Hq, gq, uq = load_qcqp_from_json(args.input)
    except Exception as e:
        print(f"Error loading QCQP problem: {e}", file=sys.stderr)
        sys.exit(1)

    # Solve
    try:
        solution = solve_qcqp(H, g, A, b, C, d, lb, ub, Hq, gq, uq,
                              solver=args.solver, verbose=args.verbose)
    except Exception as e:
        print(f"Error solving QCQP: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

    # Save
    try:
        save_solution_to_json(solution, args.output)
        print(f"Solution saved to {args.output} (solver: {solution['solver_used']}, "
              f"status: {solution['status']}, time: {solution['solve_time']:.6f}s)")
    except Exception as e:
        print(f"Error saving solution: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
