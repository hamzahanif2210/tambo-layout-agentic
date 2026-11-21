# Agentic AI loop search for TAMBO layout optimization (proof of concept)

An optimizer (differential evolution, local search, or an LLM) writes a short *plan*
(units per zone, spacing, orientation, pattern). A deterministic placer turns the plan into
unit positions on the terrain, the layout is re-scored exactly on existing CORSIKA surface
particles, and a gated reward is fed back.

## Layout
- `scripts/` Julia/Python: placer (`Placer.jl`), objective (`Objective.jl`), DE driver (`12_de_driver.jl`),
  local search (`16_local_search.jl`), event bootstrap (`14_bootstrap.jl`), LLM driver (`17_llm_driver.py`), figures.
- `slurm/` sbatch wrappers (serial_requeue, single core).
- `analysis/` bootstrap tables and notes. `plans/` example plans. `talk/` slides (LaTeX + PDF).

## Notes
- Needs the TAMBO simulation packages (TamboSim, LlhamaFit.jl) and the official CORSIKA
  simulations (read-only); paths are cluster specific.
- The LLM driver reads an API key from `~/.nrp_api` (not included). Endpoint: NRP OpenAI-compatible API.
- Data splits: train / select / test; the test split is unused so far.
- Site: the `volcano_valley` simulation dataset.
