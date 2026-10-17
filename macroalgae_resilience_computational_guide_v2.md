# Computational Guide: Simulating Resilience and Domains of Attraction in a Coral–Macroalgae Community

### Version 2 — with bifurcation testing, early-noise diagnostics, a less conservative Lyapunov estimate, and identifiability tooling

> **How to read this document.** Every technical section keeps its original math and procedure. Wherever this version adds or changes something, you'll see a box like this:
>
> **Plain-language: why this step exists.** A short explanation of the idea without notation, aimed at someone who wants to understand *why* a step is there before deciding whether to build it.

---

## 1. Project idea

This project creates a synthetic marine-ecology data set and uses it to study whether a disturbed coral–macroalgae community returns to a desired state or shifts toward macroalgal dominance.

The project is educational and computational. It does **not** claim to reproduce a particular reef. Instead, it gives researchers a transparent workflow for testing model structure, data requirements, parameter estimation, Lyapunov functions, domains of attraction, disturbances, and uncertainty before using field data.

The central question is:

> After a disturbance, which combinations of coral cover, macroalgae, herbivory, and nutrients lead to recovery toward a coral-dominated state?

**What's new in this version.** The original workflow treated four things as settled that are actually open questions or missed opportunities:

1. It assumed two stable "valleys" (coral-dominated and algae-dominated) exist, without testing whether the model actually produces them.
2. Its Lyapunov "safe region" estimate used the simplest possible shape, which is almost certainly far smaller than the true safe region — and the guide didn't say by how much.
3. It treated randomness (storms varying in strength, noisy weather, natural population fluctuation) as an optional add-on at the very end, even though resilience is fundamentally about how a system behaves *under* disturbance and noise.
4. It warned that not all parameters could be pinned down from data, but didn't give a way to find out *which ones*.

This version adds four corresponding pieces: a **bifurcation analysis** phase, an improved **Lyapunov-equation-based** candidate with a reported conservatism gap, an **early stochastic pass**, and **identifiability diagnostics** (profile likelihood and Fisher information). Everything else from the original guide is preserved.

---

## 2. Conceptual model

Imagine the reef as a marble moving across a landscape. A coral-dominated community and a macroalgae-dominated community may each be valleys where the marble can settle. A disturbance pushes the marble to a new location.

- If the marble remains in the coral valley, the reef recovers toward coral dominance.
- If it crosses the boundary between valleys, the reef may shift toward macroalgal dominance.
- The domain of attraction is the set of starting points that end in one particular valley.
- A Lyapunov function is a score that should decrease as the system moves toward a chosen target valley.

> **Plain-language: the marble metaphor has a hidden assumption.** The two-valley picture assumes the landscape actually *has* two valleys. It's equally possible the landscape is just one long slope that gets steeper as stress increases — meaning the reef degrades gradually rather than "tipping." Both pictures are debated among reef scientists. Section 8.5 below is what lets your model *tell you* which picture it produces, instead of you assuming it going in.

This is a simplified representation of ecological resilience. Real reefs are spatial, seasonal, stochastic, and subject to model uncertainty. Those complexities can be added after the basic model works.

---

## 3. Research questions

1. Does the model contain one or several ecologically meaningful equilibria?
2. Is the coral-dominated equilibrium locally stable?
3. Is the macroalgae-dominated equilibrium also stable?
4. **(New)** Does the number and stability of equilibria change as key parameters (nutrients, grazing, fishing) vary — i.e., does the model exhibit a genuine bifurcation, or just a smooth gradient of outcomes?
5. Which initial conditions lead to coral recovery?
6. Which initial conditions lead to macroalgal dominance?
7. How does a storm, heatwave, nutrient pulse, or herbivore loss change the outcome?
8. Does a candidate Lyapunov function decrease inside a recovery region, and how large is that region relative to the true (empirical) basin?
9. **(New)** How does adding realistic noise change recovery probability and the effective size of the safe region, especially near the basin boundary?
10. How much does the estimated recovery region change when parameters are uncertain?
11. Which ecological processes most strongly control resilience?
12. **(New)** Which parameters (or combinations of parameters) can actually be estimated from realistic monitoring data, and which cannot?
13. Which management intervention most increases the recovery region?

---

## 4. State variables

Begin with four state variables:

\[
x(t) = (C(t), M(t), H(t), N(t))
\]

| Symbol | Meaning | Suggested unit |
|---|---|---|
| \(C\) | Live coral cover | proportion from 0 to 1 |
| \(M\) | Macroalgal cover | proportion from 0 to 1 |
| \(H\) | Herbivore biomass or grazing capacity | standardized biomass |
| \(N\) | Nutrient concentration | standardized concentration |

The first version can use only \(C\), \(M\), and \(H\). Nutrients can initially be treated as an external forcing variable and later become a dynamic state variable.

Because coral and macroalgae compete for space, impose the ecological constraint:

\[
C + M \leq 1.
\]

If turf algae, rubble, and bare substrate are included, define available space as:

\[
B = 1 - C - M - T - R,
\]

where \(T\) is turf cover and \(R\) is rubble or other substrate. For a beginner project, retain only coral, macroalgae, and an implicit unoccupied-space term.

---

## 5. Mechanistic dynamical model

Use the general form:

\[
\frac{dx}{dt} = f(x, \theta, u(t)),
\]

### 5.1 Coral equation

\[
\frac{dC}{dt} = r_C C\left(1 - \frac{C}{K_C}\right)
- a_{CM}CM
- \mu_C C
+ r_{rec}B.
\]

### 5.2 Macroalgae equation

\[
\frac{dM}{dt} = r_M M\left(1 - \frac{M}{K_M}\right)
+ \beta_N N M
- gHM
- a_{MC}CM
- \mu_M M.
\]

A saturating grazing term can replace the linear \(gHM\) term:

\[
G(M,H) = G_{max}H\frac{M}{K_G+M}.
\]

### 5.3 Herbivore equation

\[
\frac{dH}{dt} = r_H H\left(1 - \frac{H}{K_H}\right)
+ q_M M H
- \mu_H H
- F_H H.
\]

### 5.4 Nutrient equation

\[
\frac{dN}{dt} = I_N - \lambda_N N - u_N(t).
\]

For the first computational experiment, hold \(N\) constant in low-, medium-, and high-nutrient scenarios.

---

## 6. Parameter table

| Parameter | Meaning | Example starting value | Sensitivity priority |
|---|---|---:|---|
| \(r_C\) | Coral growth rate | 0.20 | High |
| \(K_C\) | Coral capacity | 0.70 | Medium |
| \(r_M\) | Macroalgal growth rate | 0.45 | High |
| \(K_M\) | Macroalgal capacity | 0.80 | Medium |
| \(a_{CM}\) | Macroalgal effect on coral | 0.35 | High |
| \(a_{MC}\) | Coral effect on macroalgae | 0.10 | Medium |
| \(\mu_C\) | Coral mortality | 0.05 | High |
| \(\mu_M\) | Macroalgal mortality | 0.12 | Medium |
| \(g\) | Grazing coefficient | 0.40 | High |
| \(r_H\) | Herbivore growth rate | 0.15 | Medium |
| \(K_H\) | Herbivore capacity | 1.00 | Medium |
| \(\mu_H\) | Herbivore mortality | 0.08 | High |
| \(F_H\) | Fishing pressure | 0.00–0.20 | High |
| \(\beta_N\) | Nutrient-growth effect | 0.20 | High |
| \(N\) | Nutrient scenario | 0.5–2.0 | High |
| \(r_{rec}\) | Coral recruitment | 0.05 | High |

These are starting values chosen to produce interpretable simulations, not field estimates.

---

## 7. Ecological constraints

\[
C \geq 0, \quad M \geq 0, \quad H \geq 0, \quad N \geq 0, \qquad C + M \leq 1.
\]

Implementation strategies:

1. Small time step, clip invalid values after every step.
2. Transformed (e.g., logistic) coordinates.
3. Explicit space limitation in the equations.
4. Reject parameter sets that create biologically impossible trajectories.

---

## 8. Finding equilibria

An equilibrium \(x^*\) satisfies \(f(x^*) = 0\).

### Numerical procedure

1. Create many possible starting guesses.
2. Use a nonlinear root solver from each guess.
3. Keep solutions that are nonnegative and ecologically feasible.
4. Remove duplicate solutions.
5. Label each equilibrium according to its composition.

Do not assume that the model has alternative stable states. Check whether they actually arise for the selected parameters — Section 8.5 turns this check into a full analysis.

---

## 8.5 Bifurcation and regime structure (NEW — insert as Phase 1.5)

This is the direct fix for the "assumed bistability" weakness. Instead of checking equilibria at one fixed parameter set, sweep a control parameter and track how the equilibria change.

### Procedure

1. Choose one control parameter at a time — good first choices are grazing \(g\), nutrient level \(N\), or fishing pressure \(F_H\), since these are the parameters most directly tied to management levers.
2. Step the control parameter across a realistic range (e.g., \(g\) from 0 to 0.8 in small increments).
3. At each step, find all equilibria using the previous step's equilibria as the initial guesses for the root solver (this "continuation" trick keeps you on the same solution branch instead of jumping around).
4. At each equilibrium, compute the Jacobian eigenvalues (Section 9) and record whether it is stable or unstable.
5. Plot equilibrium coral cover \(C^*\) against the control parameter. This is a **bifurcation diagram**.

### What to look for

- **A single smooth curve** of stable equilibria as the parameter varies → no bistability; the reef degrades gradually as conditions worsen. This is itself an important, reportable result.
- **An S-shaped or folded curve**, where for some range of the control parameter there are two stable branches (coral-dominated and algae-dominated) separated by an unstable middle branch → this is a **saddle-node (fold) bifurcation**, the actual mathematical signature of hysteresis and tipping points.
- The parameter value where the fold appears/disappears is the **tipping threshold** for that driver — report it explicitly (e.g., "the model predicts a fold bifurcation in grazing near \(g \approx 0.18\)").

### Output

A new figure ("Figure 0"): equilibrium coral cover vs. control parameter, with solid lines for stable branches and dashed lines for unstable branches, for each of \(g\), \(N\), and \(F_H\) swept independently.

> **Plain-language: how and why.**
> **How:** you re-run the "find the resting point" calculation many times while slowly turning one dial (say, water nutrient level), and plot where the resting point ends up each time.
> **Why:** the marble-and-two-valleys story is only true if the landscape actually folds into two valleys for realistic conditions. By turning the dial and watching what happens, the model shows you its own landscape instead of you assuming one. If the plot shows one smooth curve, there's no tipping point in this model — the reef just gets steadily worse as pressure increases, and the entire "domain of attraction" framing should be reported with that caveat, not applied by default. If the plot shows a fold, you've found real evidence of tipping behavior, and you also get the exact threshold where it happens — genuinely useful for management, since it tells you the dial value to stay under.

---

## 9. Local stability using the Jacobian

The Jacobian matrix is the matrix of partial derivatives:

\[
J(x) = \frac{\partial f}{\partial x}.
\]

At an equilibrium \(x^*\), calculate the eigenvalues of \(J(x^*)\).

- All eigenvalues have negative real parts: locally asymptotically stable.
- At least one eigenvalue has a positive real part: unstable.
- An eigenvalue has zero real part: the simple test is inconclusive (this often happens *exactly at* a bifurcation point found in Section 8.5).

> **Plain-language: how and why.** The Jacobian is a "if I nudge the system slightly, does it spring back or run away?" test. You nudge each state variable a tiny bit and ask, mathematically, whether the nudge grows or shrinks over time. If every possible nudge shrinks back down, the equilibrium is a stable valley bottom. This only tells you about *very small* nudges near the equilibrium — it says nothing about what happens after a real storm, which is why you still need the basin-mapping and Lyapunov steps below.

---

## 10. Simulating trajectories

For each initial state \(x_0\), numerically integrate \(\dot{x}=f(x)\) using a standard ODE solver (explicit or adaptive Runge–Kutta; SciPy or Julia's DifferentialEquations ecosystem).

For a grid of initial conditions:

1. Choose values of \(C_0\), \(M_0\), and \(H_0\).
2. Reject states violating ecological constraints.
3. Simulate each trajectory for a long enough time.
4. Classify its final state.
5. Plot the result in a phase plane (coral on x-axis, macroalgae on y-axis; blue = coral-dominated, green = mixed, red = algae-dominated, grey = unresolved).

This empirical basin map is a useful first approximation to the domains of attraction.

---

## 10.5 Early stochastic pass (NEW — insert as Phase 3.5, right after basin mapping)

The original guide put randomness at the very end (Phase 7). This version moves a *minimal* stochastic version in immediately after the deterministic basin map, because resilience is fundamentally about behavior under disturbance and noise — a purely deterministic basin can be misleading on its own.

### Procedure

1. Convert the ODE to a simple stochastic differential equation (SDE) by adding noise to each equation, e.g. \(dC = f_C(x)\,dt + \sigma_C\,dW_C\), and similarly for \(M\) and \(H\). Start with additive Gaussian noise; multiplicative noise (noise proportional to the state) is a reasonable second step since cover proportions can't go negative.
2. Integrate using the Euler–Maruyama method (the stochastic equivalent of a simple Euler step) — no need for anything more advanced at this stage.
3. Run many replicate simulations from the same initial condition and noise level.
4. Compute two diagnostics:
   - **Early-warning signals:** as a control parameter approaches the bifurcation threshold found in Section 8.5, track the variance and lag-1 autocorrelation of the simulated trajectory over a moving window. Both are expected to rise near a fold bifurcation ("critical slowing down").
   - **Mean first-passage time / transition probability:** for a fixed noise level, estimate how long it typically takes (or how likely it is within a fixed window) for a trajectory started near the coral-dominated equilibrium to cross into the algae-dominated basin, purely due to noise — with no discrete disturbance event.

### Output

Two new figures: (a) variance/autocorrelation vs. distance from the bifurcation threshold; (b) transition probability vs. noise intensity, for a few starting points near the basin boundary identified in Section 10.

> **Plain-language: how and why.**
> **How:** instead of a clean, predictable trajectory, you now let the model take a slightly random walk at every time step — like adding realistic day-to-day wobble to fish populations or wave conditions — and see where it ends up, many times over.
> **Why:** a resting point can be "stable" in the strict mathematical sense but still be fragile in practice, if a little bit of everyday randomness is enough to nudge the system across the boundary into the other valley. This is important for real management: a reef sitting *near* the boundary of its basin might look safe in the clean, deterministic picture but actually have a meaningful chance of tipping over just from ordinary environmental noise. Early-warning signals (rising variance and "stickiness" in the data) are also one of the few things field researchers can actually measure in advance of a real tipping event, so producing them here previews what a monitoring program should watch for.

---

## 11. Defining the target state

Select a target equilibrium, e.g. \(x^* = (C^*, M^*, H^*)\), with high coral cover, low macroalgal cover, and sufficient herbivore biomass. The target must be defined from the model, not only from a visual preference — it should be an equilibrium confirmed in Sections 8 and 8.5.

---

## 12. Constructing a Lyapunov candidate

Define the error from the target: \(e = x-x^*\).

A simple candidate is:

\[
V(x) = w_C(C-C^*)^2 + w_M(M-M^*)^2 + w_H(H-H^*)^2,
\]

or in vector form, \(V(x)=(x-x^*)^T P(x-x^*)\), with \(P\) positive-definite. The original guide starts with a diagonal \(P\) (weights \(w_C, w_M, w_H\), often all set to 1).

> **Plain-language: why the diagonal choice is a weak starting point.** A diagonal \(P\) treats the three variables as if they don't interact — as if coral, algae, and herbivores each independently pull the "distance from target" score up or down. In reality they're coupled (coral loss and algae gain are linked through competition, for instance), so a diagonal shape is a crude, round-ish "safe zone" drawn around the target when the true safe zone is more like a tilted, stretched oval. Section 13.5 below replaces the guess with a shape computed directly from the model's own dynamics.

---

## 13. Computing the Lyapunov derivative

\[
\dot{V}(x) = \nabla V(x)^T f(x)
= 2w_C(C-C^*)\frac{dC}{dt} + 2w_M(M-M^*)\frac{dM}{dt} + 2w_H(H-H^*)\frac{dH}{dt}.
\]

- \(V=0\) at the target; \(V>0\) elsewhere.
- \(\dot V<0\) means the score is decreasing (moving toward the target); \(\dot V>0\) means it is increasing.

A strict global proof that \(\dot V<0\) everywhere may be impossible. The practical goal is to identify a useful region where the conditions hold.

---

## 13.5 A better Lyapunov shape, and measuring how conservative it is (NEW)

This directly addresses the "Lyapunov region is probably tiny and we don't say by how much" weakness.

### Step 1 — Solve for \(P\) instead of guessing it

1. Linearize the system at the target equilibrium: compute the Jacobian \(J = \partial f/\partial x\) at \(x^*\) (already computed in Section 9).
2. Choose a simple positive-definite matrix \(Q\) (the identity matrix is a reasonable default).
3. Solve the continuous-time Lyapunov equation for \(P\):

\[
J^T P + P J = -Q.
\]

Standard numerical solvers handle this directly (e.g., a Lyapunov-equation solver in a scientific computing library). This is not a guess — it is the shape guaranteed to make \(\dot V < 0\) in some neighborhood of the target, given how the model actually behaves near that point, and it will generally produce a noticeably larger and better-shaped sublevel set than an arbitrary diagonal \(P\).

### Step 2 — Search for the sublevel set as before

Use \(\Omega_c = \{x: V(x) \leq c\}\) and search for the largest \(c\) such that \(V>0\) away from the target, \(\dot V<0\) inside \(\Omega_c\), the set stays within the feasible state space, trajectories don't leave it, and no unwanted equilibrium lies inside it (unchanged from the original guide's Section 14).

### Step 3 — Report the conservatism gap explicitly

1. Compute the area (or volume, in 3 state variables) of the empirical basin from Section 10.
2. Compute the area/volume of \(\Omega_{c^*}\).
3. Report the ratio as a headline number: *"The certified Lyapunov region covers approximately X% of the empirically estimated basin."*
4. Repeat with both the diagonal \(P\) and the Lyapunov-equation \(P\), so the improvement from Step 1 is itself quantified.

### Optional extension

Sum-of-squares (SOS) programming can search over higher-order polynomial Lyapunov candidates (not just quadratics), typically certifying a substantially larger region than any quadratic form. This is more advanced tooling (e.g., `SumOfSquares.jl`, or SOS-capable convex solvers in Python) and is worth trying only after Steps 1–3 above are working, since it replaces "solve one linear equation" with "solve a convex optimization problem."

> **Plain-language: how and why.**
> **How:** instead of drawing an arbitrary round "safe zone" around the target and hoping it works, you ask the model itself — via a standard piece of linear algebra called the Lyapunov equation — for the best-fitting oval shape given how the system actually behaves near the target. You then find the largest version of that shape for which you can still guarantee the score keeps decreasing everywhere inside it, and you directly compare its size to the messier, simulation-based safe zone from Section 10.
> **Why:** the original approach's biggest silent problem was that its "provably safe" region was probably much smaller than the real one, without saying so — which could wrongly suggest a management intervention isn't helping when it actually is, just outside the tiny certified zone. Solving for the matrix instead of guessing it usually closes much of that gap for free, and reporting the ratio turns a hidden weakness into an honest, quantified result: "our certified guarantee covers this fraction of what the simulations suggest is actually safe."

---

## 14. Estimating a Lyapunov sublevel set

(Unchanged from the original — see Section 13.5 above for the updated candidate and the new conservatism-ratio reporting step.)

\[
\Omega_c \subseteq \text{true domain of attraction}
\]

may be guaranteed, while the reverse inclusion usually is not known.

---

## 15. Empirical domain-of-attraction mapping

For each initial state \(x_0\): simulate the model, calculate the final distance to each equilibrium, assign the trajectory to the nearest attractor if convergence is clear, mark it unresolved if it remains near a boundary or does not converge, and repeat on a dense grid or with random initial states.

Compare the empirical basin with the Lyapunov sublevel set (now with the quantified ratio from Section 13.5).

---

## 16. Simulating disturbances

### 16.1 Instantaneous shock

\[
C(t_d^+)=(1-s_C)C(t_d^-), \qquad M(t_d^+) = M(t_d^-)+s_M, \qquad H(t_d^+)=(1-s_H)H(t_d^-).
\]

### 16.2 Time-varying forcing

\[
N(t)=N_0 + A_N\exp\left[-\frac{(t-t_d)^2}{2\sigma_N^2}\right].
\]

### Disturbance experiments

| Experiment | Coral loss | Herbivore loss | Nutrient pulse |
|---|---:|---:|---:|
| Mild storm | Low | Low | None |
| Strong storm | High | Medium | None |
| Heatwave | Medium | Low | None |
| Pollution event | Low | Low | High |
| Fishing pressure | None | High | None |
| Compound event | High | High | High |

For every scenario, record whether the trajectory returns to the target, approaches another attractor, or remains unresolved — and, using the noise model from Section 10.5, also record the *probability* of recovery under repeated stochastic replicates of the same disturbance, not just a single deterministic outcome.

---

## 17. Synthetic observations

Generate observed data using \(y_i = x(t_i) + \varepsilon_i\), with distributions tailored to variable type (beta for cover proportions, lognormal for biomass, normal for temperature, detection error for counts). Add missing observations and irregular survey intervals, but keep a clean version for validation.

---

## 18. Synthetic sampling design

Create several virtual reef sites spanning different baseline conditions (high coral/low algae, mixed, low coral/high algae, high/low nutrients, protected/fished). For each: simulate a burn-in period, sample monthly or quarterly, apply random disturbances, continue sampling through recovery, and record both true and observed values, with multiple replicates per scenario.

---

## 19. Parameter estimation exercise

Treat synthetic observations as if they were field data. Estimate parameters by minimizing:

\[
\mathcal{L}(\theta)=\sum_i \left\|y_i-x(t_i;\theta)\right\|^2,
\]

subject to ecological constraints, using nonlinear least squares, maximum likelihood, Bayesian inference, approximate Bayesian computation, or particle filtering.

Then compare the true parameters, estimated parameters, model predictions, and the known true trajectory.

---

## 19.5 Identifiability diagnostics (NEW)

This is the direct fix for "we warn some parameters won't be identifiable, but don't say how to find out which."

### Step 1 — Profile likelihood

For each parameter \(\theta_j\):

1. Fix \(\theta_j\) at a sequence of values spanning a plausible range around its fitted value.
2. At each fixed value, re-optimize every *other* parameter to minimize \(\mathcal{L}(\theta)\).
3. Plot the resulting best-possible \(\mathcal{L}\) against the fixed value of \(\theta_j\) — this is the **profile likelihood** for that parameter.
4. A sharp, well-defined minimum means \(\theta_j\) is practically identifiable from this data. A flat or very wide profile means the data barely constrains it, even if the optimizer reports a single "best" value.

### Step 2 — Fisher Information Matrix (FIM)

1. At the fitted parameters, numerically compute the Hessian (second-derivative matrix) of \(\mathcal{L}(\theta)\) via finite differences.
2. Compute its eigenvalues and eigenvectors.
3. Large eigenvalues correspond to "stiff" directions in parameter space (well constrained by the data); very small eigenvalues correspond to "sloppy" directions (poorly constrained). Ecological ODE models very commonly show a handful of stiff directions and many sloppy ones — this is worth reporting as a named result ("sloppy model" behavior), not treated as a bug.
4. The eigenvectors of the small eigenvalues often reveal that it's not individual parameters but specific *combinations* (e.g., a ratio like \(a_{CM}/\mu_C\)) that are actually identifiable, even when the individual parameters are not.

### Output

A new table: each parameter with its profile-likelihood shape (sharp/flat) and its FIM ranking (stiff/sloppy), plus a note on which parameter *combinations*, if any, are identifiable even when individual parameters aren't.

> **Plain-language: how and why.**
> **How:** for each unknown quantity in the model, you temporarily "pin" it to a range of possible values and let the computer find the best-fitting version of everything else at each pinned value. If changing the pinned value barely changes how well the model fits the data, that quantity isn't really being measured by your data — the optimizer just picked *a* number, not necessarily the *right* number. The Fisher Information step does something similar but faster, using calculus to directly measure how sharply the fit responds to changes in each parameter, and in which directions in parameter-space (which might mix several parameters together) the response is sharp versus flat.
> **Why:** an optimizer will always hand you a number for every parameter, whether or not the data actually contains enough information to justify it — it has no way of telling you "I'm not really sure about this one." Without this step, a fitted parameter table looks equally confident everywhere, which is misleading. With it, you can honestly say things like "grazing rate and coral mortality are well constrained by this monitoring design, but macroalgal mortality is not — more frequent sampling or better replication would be needed before trusting that value," which is exactly the kind of guidance a real monitoring program needs before it invests money in field data collection.

---

## 20. Validation plan

### Numerical validation
- Smaller solver time steps give similar trajectories.
- Mass or space constraints hold.
- Equilibrium residuals are near zero.
- Repeated runs with the same seed match.

### Dynamical validation
- Simulated equilibria match known analytical cases.
- Jacobian stability result confirmed by small perturbation simulations.
- Bifurcation diagram (Section 8.5) is stable under finer parameter-step resolution.
- Trajectories classified as convergent actually converge.

### Statistical validation
- Hold out some synthetic observations.
- Fit the model to training data; test on held-out data.
- Compare forecast errors across scenarios.
- Cross-check profile-likelihood conclusions (Section 19.5) against a second estimation method (e.g., Bayesian posterior width).

### Ecological validation
- Signs of effects are ecologically reasonable.
- Growth, mortality, and grazing rates compared with literature ranges.
- Discuss which processes the model omits.

---

## 21. Uncertainty analysis

### Parameter uncertainty
\[
\theta^{(k)} \sim p(\theta).
\]

### Initial-state uncertainty
\[
x_0^{(k)} \sim p(x_0).
\]

### Environmental uncertainty
Vary temperature, nutrients, disturbances, and other external drivers — this now also includes the noise-intensity parameter from Section 10.5 as an uncertain quantity in its own right.

For every sampled scenario, calculate: final community state, time to recovery, maximum macroalgal cover, probability of coral recovery, probability of transition to macroalgal dominance, estimated Lyapunov minimum or maximum, and whether the trajectory crosses the empirical basin boundary.

Report distributions, not only one prediction.

---

## 22. Sensitivity analysis

Useful outputs: probability of coral recovery, size of the estimated recovery region, time to return to target, maximum macroalgal cover, coral loss after disturbance, and (new) the location of the bifurcation threshold from Section 8.5.

Methods: one-at-a-time perturbation, local derivatives, Latin hypercube sampling, Sobol sensitivity indices, Morris screening, partial rank correlation.

The simulation should determine which parameters matter most rather than assume it.

---

## 23. Management experiments

### Herbivore protection
Reduce \(F_H\).

### Nutrient reduction
Reduce \(I_N\) or \(N\).

### Coral restoration
Increase \(C\) or \(r_{rec}\) after a disturbance.

### Macroalgal removal
Reduce \(M\) at a selected intervention time.

### Combined intervention
Apply multiple actions together.

For each intervention, estimate recovery probability, time to recovery, macroalgal peak, required intervention intensity, change in the estimated domain of attraction, change in distance from the bifurcation threshold (new — how much "margin" the intervention buys relative to Section 8.5's tipping point), and cost/effort if available.

> Does the intervention move the system back inside a recovery region — and how far from the tipping threshold does it move it?

---

## 24. Suggested project structure

```text
macroalgae-resilience/
├── README.md
├── environment.yml or pyproject.toml
├── data/
│   ├── synthetic_clean/
│   ├── synthetic_observed/
│   └── metadata/
├── src/
│   ├── model.py
│   ├── parameters.py
│   ├── simulate.py
│   ├── bifurcation.py          # NEW — Section 8.5
│   ├── stochastic.py           # NEW — Section 10.5
│   ├── disturbances.py
│   ├── equilibria.py
│   ├── stability.py
│   ├── lyapunov.py             # updated — solves Lyapunov equation, Section 13.5
│   ├── basin_mapping.py
│   ├── synthetic_data.py
│   ├── estimation.py
│   ├── identifiability.py      # NEW — Section 19.5
│   └── uncertainty.py
├── notebooks/
│   ├── 01_model_baseline.ipynb
│   ├── 02_equilibria_and_stability.ipynb
│   ├── 02b_bifurcation_analysis.ipynb      # NEW
│   ├── 03_domain_of_attraction.ipynb
│   ├── 03b_stochastic_basin.ipynb          # NEW
│   ├── 04_lyapunov_analysis.ipynb
│   ├── 05_disturbance_experiments.ipynb
│   ├── 06_parameter_estimation.ipynb
│   ├── 06b_identifiability.ipynb           # NEW
│   └── 07_uncertainty_and_sensitivity.ipynb
├── tests/
│   ├── test_model.py
│   ├── test_constraints.py
│   └── test_lyapunov.py
└── figures/
```

---

## 25. Recommended implementation order

### Phase 1: Deterministic model
1. Implement the coral–macroalgae–herbivore equations.
2. Simulate one trajectory.
3. Plot each state through time.
4. Check biological constraints.

### Phase 1.5: Bifurcation and regime structure *(NEW)*
1. Sweep a control parameter (\(g\), \(N\), or \(F_H\)).
2. Track equilibria and their stability across the sweep.
3. Build the bifurcation diagram; identify fold points if present.
4. Report whether the model is genuinely bistable or only smoothly graded.

### Phase 2: Equilibria and stability
1. Find equilibria from many initial guesses.
2. Calculate Jacobian eigenvalues.
3. Perturb each equilibrium numerically.
4. Confirm local stability through simulations.

### Phase 3: Basin mapping
1. Generate a grid of initial states.
2. Simulate every grid point.
3. Classify long-term outcomes.
4. Plot empirical domains of attraction.

### Phase 3.5: Early stochastic pass *(NEW — moved up from the end)*
1. Add noise to the model; integrate with Euler–Maruyama.
2. Compute early-warning signals near the bifurcation threshold.
3. Estimate transition probabilities near the basin boundary.

### Phase 4: Lyapunov analysis
1. Choose a target equilibrium.
2. Solve the Lyapunov equation for \(P\) at the linearization (not just a diagonal guess).
3. Calculate \(V\) and \(\dot V\).
4. Search for a valid sublevel set.
5. Compare it with the simulation-based basin and report the conservatism ratio.

### Phase 5: Disturbances
1. Add instantaneous shocks.
2. Add heatwave and nutrient forcing.
3. Simulate single and compound events, deterministically and stochastically.
4. Measure recovery probability and time.

### Phase 6: Synthetic observations
1. Add measurement error.
2. Add missing observations.
3. Fit the model to the observed data.
4. Compare estimated and true parameters.

### Phase 6.5: Identifiability diagnostics *(NEW)*
1. Compute profile likelihoods for each parameter.
2. Compute the Fisher Information Matrix and its eigenstructure.
3. Report which parameters (or combinations) are identifiable from this sampling design.

### Phase 7: Modern extensions
1. Add richer environmental stochasticity (beyond the minimal pass in Phase 3.5).
2. Try sum-of-squares Lyapunov candidates.
3. Try data-driven or hybrid modelling.
4. Add spatial patches.
5. Compare interventions using the bifurcation margin (Section 23).

---

## 26. Core outputs

### Figures
1. Coral, macroalgae, and herbivore trajectories.
2. Phase-plane trajectories.
3. Equilibrium locations.
4. **Bifurcation diagram** (equilibrium vs. control parameter, stable/unstable branches). *(NEW)*
5. Empirical domains of attraction.
6. Lyapunov contours \(V(x)=c\), for both the diagonal and Lyapunov-equation \(P\). *(updated)*
7. Regions where \(\dot V<0\), with the empirical-vs-certified area/volume ratio annotated. *(updated)*
8. **Early-warning signals** (variance/autocorrelation) vs. distance from the bifurcation threshold. *(NEW)*
9. **Stochastic transition probability** vs. noise intensity near the basin boundary. *(NEW)*
10. Disturbance and recovery trajectories.
11. Recovery probability under disturbance size (deterministic and stochastic).
12. Sensitivity rankings.
13. **Profile likelihood curves** for each parameter. *(NEW)*
14. Intervention comparison, including margin from the bifurcation threshold.

### Tables
1. Parameter definitions and values.
2. Equilibrium states and eigenvalues.
3. Bifurcation thresholds by control parameter. *(NEW)*
4. Disturbance scenarios.
5. Recovery probabilities (deterministic and stochastic).
6. Parameter-identifiability results (profile likelihood + FIM ranking + identifiable combinations). *(updated)*
7. Sensitivity results.
8. Lyapunov conservatism ratio (diagonal vs. Lyapunov-equation \(P\)). *(NEW)*
9. Model limitations.

---

## 27. Interpretation rules

Prefer:

> Under the specified model and parameter ranges, initial states inside the estimated region were predicted to recover toward the coral-dominated equilibrium.

Avoid:

> The reef will recover whenever it is inside this exact boundary.

Prefer:

> The Lyapunov set is a certified or numerically supported inner estimate, covering approximately X% of the empirically mapped basin.

Avoid:

> The Lyapunov set is the true ecological basin.

Prefer:

> A bifurcation analysis found a fold in the grazing parameter near \(g \approx g_{crit}\); below this value the model no longer exhibits alternative stable states.

Avoid:

> This model proves reefs have tipping points. *(It shows this particular model, under these parameters, can — not that real reefs necessarily do.)*

Prefer:

> Profile-likelihood analysis indicates that grazing rate and coral mortality are well constrained by this sampling design, while macroalgal mortality is not.

Avoid:

> The fitted parameters are the true parameters. *(A single best-fit number says nothing about how confidently it was estimated.)*

Prefer:

> The simulation suggests that herbivory and nutrient enrichment strongly influence the basin boundary.

Avoid:

> Herbivory alone controls resilience.

The model is a structured hypothesis, not a replacement for field validation.

---

## 28. Limitations and extensions

### Limitations of the basic model
- Assumes a well-mixed reef.
- May omit seasonal cycles.
- May treat parameters as constant.
- May not represent species-specific interactions.
- May ignore spatial refuges for macroalgae.
- May ignore larval connectivity.
- May not represent abrupt tipping behavior accurately *(now directly testable via Section 8.5, rather than assumed)*.
- May use a simplified disturbance model.
- A deterministic Lyapunov function does not fully represent stochastic ecosystems *(now partially addressed by Section 10.5, though only with a minimal noise model)*.
- Profile-likelihood and FIM diagnostics assume the model structure is correct; they diagnose *data* limitations, not model-misspecification.

### Extensions
- Add spatial reef patches.
- Add explicit turf algae, rubble, and substrate.
- Add temperature-dependent growth and mortality.
- Use richer stochastic differential equations (state-dependent noise, colored noise).
- Include multiple coral functional groups.
- Model fish movement and fishing spatially.
- Use remote-sensing observations.
- Estimate parameters with Bayesian methods (and compare posterior width directly against the profile-likelihood results).
- Use reachability analysis for uncertainty bounds.
- Learn a Lyapunov function using sum-of-squares or neural methods, and compare its certified region against the linear-Lyapunov-equation result.
- Compare the model against real monitoring data.

---

## 29. Minimal pseudocode

```python
# 1. Define parameters and ecological equations
params = make_parameters()

def rhs(t, x, params, forcing):
    C, M, H, N = x
    dC = coral_growth(C, M, H, N, params, forcing)
    dM = macroalgae_growth(C, M, H, N, params, forcing)
    dH = herbivore_dynamics(C, M, H, N, params, forcing)
    dN = nutrient_dynamics(C, M, H, N, params, forcing)
    return [dC, dM, dH, dN]

# 2. Simulate baseline and disturbed trajectories
x0 = [0.60, 0.15, 0.70, 0.50]
trajectory = integrate(rhs, x0, params)
disturbed = apply_disturbance(trajectory, event="storm")

# 3. Find equilibria and test local stability
equilibria = find_equilibria(rhs, params)
for x_star in equilibria:
    J = numerical_jacobian(rhs, x_star, params)
    eigenvalues = eigvals(J)

# 3.5 NEW: bifurcation sweep over a control parameter
bifurcation_diagram = sweep_control_parameter(
    rhs, param_name="g", param_range=(0.0, 0.8), n_steps=200
)
fold_points = detect_fold_bifurcations(bifurcation_diagram)

# 4. Construct a Lyapunov candidate
x_target = choose_coral_dominated_equilibrium(equilibria)
J_target = numerical_jacobian(rhs, x_target, params)

# NEW: solve the Lyapunov equation instead of guessing P
Q = identity_matrix(3)
P = solve_lyapunov_equation(J_target, Q)   # J^T P + P J = -Q

def V(x):
    error = x[:3] - x_target[:3]
    return error.T @ P @ error

def Vdot(x):
    gradient = numerical_gradient(V, x)
    return gradient @ rhs(0, x, params, forcing=None)

# 5. Map the empirical basin
basin = classify_initial_conditions(rhs, equilibria, params)

# 5.5 NEW: minimal stochastic pass
sde_trajectories = simulate_sde(rhs, x0, params, noise_sigma=0.02, n_replicates=500)
early_warning = compute_variance_and_autocorrelation(sde_trajectories)
transition_prob = estimate_transition_probability(sde_trajectories, basin)

# 6. Search for a Lyapunov sublevel set
c_values = candidate_levels(V)
valid_levels = [c for c in c_values if check_Vdot_negative(c, Vdot)]
c_star = max(valid_levels)
conservatism_ratio = compute_area(V, c_star) / compute_area(basin)  # NEW

# 7. Add observation noise and fit the model
observed = generate_observations(trajectory, measurement_error=True)
fit = estimate_parameters(observed, rhs)

# 7.5 NEW: identifiability diagnostics
profiles = compute_profile_likelihoods(observed, rhs, fit)
fisher_info = compute_fisher_information_matrix(observed, rhs, fit)
identifiable_directions = eigendecompose(fisher_info)

# 8. Run uncertainty and management scenarios
results = run_parameter_ensemble(rhs, fit, disturbances, interventions)
```

---

## 30. Final deliverable

A complete computational project should include:

- a clearly documented ecological model;
- a reproducible parameter file;
- clean simulated trajectories;
- **a bifurcation analysis establishing whether the model is genuinely bistable, and under what conditions** *(new)*;
- noisy synthetic monitoring data;
- equilibrium and stability analysis;
- empirical basin maps;
- **a stochastic pass with early-warning signals and transition probabilities near the basin boundary** *(new)*;
- Lyapunov-function calculations, **including a quantified comparison between the certified region and the empirical basin** *(new)*;
- disturbance experiments;
- uncertainty and sensitivity results;
- **identifiability diagnostics (profile likelihood and Fisher information) for every estimated parameter** *(new)*;
- model-validation tests;
- figures and tables;
- a limitations section;
- instructions for replacing synthetic data with field data.

The project's main scientific contribution is methodological: it demonstrates how researchers can move from ecological assumptions to dynamic equations, from equations to equilibria, from equilibria to a *tested* (not assumed) regime structure, from that structure to recovery regions with an honestly quantified degree of certainty, and from recovery regions to disturbance, noise, and management analysis.

The final interpretation should remain cautious:

> A domain of attraction is not a permanent promise of recovery. It is a model-based estimate of which ecological states are likely to return to a target under specified assumptions, parameters, disturbances, and uncertainties — and that estimate is only as strong as the evidence, from the model's own bifurcation structure and from the data's own identifiability, that supports it.
