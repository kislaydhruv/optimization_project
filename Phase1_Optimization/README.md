# Phase 1 Optimization — Bounding Phase + Bisection Method

ME 609 — Single Variable Optimization
Files: `optimization_phase1.cpp`, `input.txt`

## How to run

```bash
g++ -O2 -std=c++17 -o phase1 optimization_phase1.cpp
./phase1
Please assign the question number
1                # <- type 1-6 here, program reads a/b/delta/eps/x0/h for that
                 #    question from input.txt automatically
```

To generate the "run 10 times with different initial guesses" table required
in the slides: open `input.txt`, edit the `x0` column for the question you're
running, save, and re-run `./phase1` with that question number. Repeat for
10 different `x0` values per question, and again after changing `delta`/`h`
to see the effect on accuracy / nEval as instructions (c) asks.

## Program structure (matches the required structure)

- `main()` — reads `input.txt`, calls `bounding_phase()`, then
  `bisection_method()`, prints/saves results.
- `bounding_phase()` — bracketing method.
- `bisection_method()` — region-elimination / gradient-based method.
- `fun_eval()` — objective_function(), switch-cased on the global `QUE`.
  Maximization problems (Q1, Q2, Q3, Q6) are converted to minimization by
  negating `f(x)` inside the switch, as required — no separate
  commenting/uncommenting needed, just answer the question-number prompt.
- `fun_der()` — other_functions(): central-difference numerical derivative,
  `(f(x+h) - f(x-h)) / 2h`.

`f(x*)` printed for Q1/Q2/Q3/Q6 is the value of the **minimized surrogate**
(`-f_original`); the true maximum is `-f(x*)`.

## Input file format

One row per question, whitespace separated:

```
QUE   a      b        delta   eps      x0     h
1     -10    0        0.5     1e-4     -5     1e-4
```

- `a, b` — interval given in the problem statement.
- `delta` — bounding-phase step size.
- `eps` — bisection tolerance on `|f'(x)|`.
- `x0` — initial guess (edit this per run).
- `h` — finite-difference step for the numerical derivative.

## Bugs found in the draft and what was fixed

Your original `draft.cpp` had the right overall approach and the function
signatures/names below were kept exactly as you wrote them. These are the
concrete mistakes found while checking the logic against the lecture slides,
and the fixes made:

1. **Q2 objective was never negated.** Q2 asks to *maximize*
   `8 + x^3 - 2x - 2e^x`, but the switch case had no `(-1)` (unlike Q1, Q3,
   Q6 which correctly negate). Fixed by adding the missing `*(-1)`.

2. **`bounding_phase()` while-loop condition was inverted**, so it never
   actually executed the expansion loop, and the point for the "rejected"
   search direction was never repositioned (only its function value was
   copied) — so `x0` could be moved to the wrong location the moment
   Δ turned out negative. Rewritten to follow the slides exactly:
   `x(k+1) = x(k) + 2^k*Δ`, keep expanding while `f(k+1) < f(k)`, and
   return the bracket `(x(k-1), x(k+1))` when it stops.

3. **`bisection_method()` took `total_func_eval` by value**, not by
   reference, so the running evaluation count from `bounding_phase()`
   never carried over and never made it back to `main()`. Changed to
   `int&` like every other function in the file.

4. **Wrong parentheses in the bisection stopping condition**:
   `abs(fun_der(z,d,total_func_eval) > e)` computes `abs()` of a boolean
   (always 0 or 1), not `|f'(z)| > e`. Fixed to
   `fabs(fun_der(...)) > e`.

5. **Repeated derivative evaluations (the "repetition" you flagged).**
   The bisection loop called `fun_der(z, ...)` up to three times per
   iteration (once in the while-condition, again inside the `if`, again
   inside the `else if`), tripling the reported `nEval` for no reason.
   Fixed to compute `der_z` once per iteration and reuse it — matches the
   slide's note that bisection "requires two function evaluations per
   iteration" (one `fun_der` call = one central difference = 2 evals).

6. **`a`, `b` in `bounding_phase()` were declared but unused.** They're
   now used to clamp the exponential search inside the interval given for
   the problem, which also protects Q4's `e^(0.5x^2)` term from overflowing
   if a bad initial guess sent the search running away.

7. **No bracket-validity check in `bisection_method()`.** Added a warning
   (not a hard failure) if `f'(a) < 0` / `f'(b) > 0` isn't satisfied, since
   that's the slide's precondition for the method.

Everything else (variable names, function names, `fun_eval`/`fun_der`
signatures, the overall two-phase structure) was left as you wrote it.

## Sample output (actual runs, all 6 questions)

Verified against a 2,000,001-point brute-force grid search over each
interval — all 6 match to 5+ decimal places.

```
=== QUE=1 : Maximize (2x-5)^4-(x^2-1)^3 in (-10,0) ===
Bracket from bounding phase : [-6.500000, -5.000000]  (nEval so far = 5)
x(0)   = -5.000000
x*     = -5.775068
f(x*)  = -41165.614407      (true max = 41165.614407)
nEval  = 62

=== QUE=2 : Maximize 8+x^3-2x-2e^x in (-2,1) ===
Bracket from bounding phase : [-1.100000, -0.900000]  (nEval so far = 3)
x(0)   = -1.000000
x*     = -0.960156
f(x*)  = -8.269478          (true max = 8.269478)
nEval  = 24

=== QUE=3 : Maximize 4x*sin(x) in (0.5,pi) ===
Bracket from bounding phase : [1.800000, 3.000000]  (nEval so far = 7)
x(0)   = 1.500000
x*     = 2.028754
f(x*)  = -7.278823          (true max = 7.278823)
nEval  = 44

=== QUE=4 : Minimize 2(x-3)^2+e^(0.5x^2) in (-2,3) ===
Bracket from bounding phase : [0.600000, 3.000000]  (nEval so far = 7)
x(0)   = 0.000000
x*     = 1.590720
f(x*)  = 7.515924
nEval  = 48

=== QUE=5 : Minimize x^2-10e^(0.1x) in (-6,6) ===
Bracket from bounding phase : [0.000000, 1.500000]  (nEval so far = 5)
x(0)   = 0.000000
x*     = 0.527069
f(x*)  = -10.263405
nEval  = 38

=== QUE=6 : Maximize 20sin(x)-15x^2 in (-4,4) ===
Bracket from bounding phase : [0.100000, 1.000000]  (nEval so far = 5)
x(0)   = 1.000000
x*     = 0.563568
f(x*)  = -5.919983          (true max = 5.919983)
nEval  = 38
```

