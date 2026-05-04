# LLM-Augmented Deep Learning for Drug Discovery: From Expert Knowledge to Implementation at Scale

**Proposal Type**: Strategic Research Initiative  
**Date**: May 2026  
**Status**: Draft for Review

---

## Executive Summary

A proof-of-concept study at our institution demonstrated that large language models (LLMs) such as Claude Code can bridge the gap between scientific publication and working implementation in ways that were previously impractical. In this study, complex inference-time steering modifications to a state-of-the-art co-folding model (Boltz-2) — modifications described only in conceptual terms in recent preprints — were proposed, implemented, evaluated, and iteratively improved across two full experimental rounds in a matter of weeks.

The implication is broader than antibody design. We are at the beginning of a paradigm shift in which **any department with domain knowledge and access to an LLM can customize, extend, and experiment with deep learning models** in drug discovery — without requiring a dedicated team of ML engineers. This proposal argues that capitalizing on this shift is not optional: the competitive advantage it confers is significant, and the cost of staying behind is measurable.

---

## 1. The New Paradigm

For the past decade, the practical gap between "a deep learning model is published" and "we are using a customized version of it for our specific problem" has been enormous. That gap was filled by ML engineers who could translate scientific ideas from papers into code, debug implementation details that were left unstated, and iterate experiments at speed. Most drug discovery departments did not have these engineers, or had too few of them.

LLMs, specifically agentic coding systems like Claude Code, are closing this gap in real time. They can:

- **Read and interpret scientific papers**, extracting the algorithmic intent behind described methods
- **Translate ideas into working code**, even when no reference implementation exists
- **Identify failure modes and propose improvements**, incorporating feedback from experimental results and domain expert guidance
- **Iterate across experimental rounds**, treating each evaluation as structured input to the next design cycle

The result is that **tasks previously requiring weeks of specialized engineering effort can now be completed in days** — and tasks previously considered too speculative to attempt are now within reach of scientists who understand the biology but not the implementation.

This is not a future possibility. It happened.

---

## 2. Proof of Concept: Steering a Co-folding Model for Antibody–Antigen Structure Prediction

### 2.1 Context

Boltz-2 is a state-of-the-art deep learning model for predicting protein complex structures, including antibody–antigen assemblies. Its outputs — predicted binding geometries between antibodies and their targets — are directly relevant to antibody drug discovery, where understanding which part of a target (the epitope) an antibody engages, and how precisely the antibody loops are positioned, determines downstream developability and efficacy.

The challenge: Boltz-2, like most diffusion-based models, predicts structures through stochastic sampling. Without intervention, it may not produce the binding geometry most relevant to a discovery program. The scientific literature describes several approaches to guide, or *steer*, diffusion models toward desired outputs at inference time — without retraining. But these descriptions were ideas and mathematical formulations, not code.

### 2.2 What Was Done

Working in two sequential rounds, the LLM (Claude Code) was given:

1. **The Boltz-2 codebase** — tens of thousands of lines of Python
2. **Expert knowledge from domain scientists** — biological constraints, knowledge of CDR loop biology, understanding of what "good" antibody–antigen structure prediction looks like
3. **Three recent scientific publications** describing relevant inference-time steering techniques:
   - FK-Diffusion Steering (Horvitz et al., arXiv:2501.06848) — Feynman-Kac importance sampling for diffusion models
   - Boltz-sample β-scaling (Suzuki & Amagasa, 2026) — pair representation modulation in latent space
   - EmbedOpt (Li et al., arXiv:2602.05285) — robust embedding-space optimization for protein diffusion

From these inputs, the LLM **designed, implemented, and evaluated 17 distinct steering strategies** across two experimental rounds, benchmarked against 47 antibody–antigen complexes with ground-truth crystal structures. Key results included:

- **Strategy Y (Hierarchical Steering)**: the only approach across all 17 strategies to achieve a statistically significant structural improvement — reducing CDR-H3 loop error by 0.30 Å (p = 0.025). CDR-H3 is the primary antibody loop for antigen engagement.
- **Strategy A (FK Particle Resampling)**: best overall docking quality improvement (+23% DockQ over baseline), with zero quality regressions — a "safety net" property that makes it suitable for production use.
- **Strategy L+ (Safe CDR3 Scaling)**: eliminated harmful side effects identified in Round 1 while maintaining positive binding quality trends — demonstrating the value of the iterative feedback loop.
- **Strategy Q (Iterative Epitope Refinement)**: best blind epitope prediction improvement, with a +60% relative gain in median prediction quality on the hardest cases where the baseline found nothing.

### 2.3 The Iterative Nature of the Work

Round 1 identified both successes and failures. The failures were as informative as the successes: some strategies degraded confidence calibration, some were harmful to structure quality, and some revealed a fundamental bottleneck (the ~68% "failure floor" where global orientation is wrong — a problem no local steering method can fix). These lessons were fed back to the LLM as structured analysis documents, which then designed Round 2 with explicit modifications targeting the identified failure modes.

This feedback loop — **propose → implement → evaluate → analyze → improve** — proceeded with a velocity that would not have been achievable with a traditional engineering team. The time from "idea described in a paper" to "implemented, tested, and analyzed" was measured in days, not months.

### 2.4 What Made This Work

Three ingredients were essential, and all three are replicable:

1. **Domain expert knowledge**: Scientists who understood CDR biology, antibody structure–activity relationships, and what metrics mattered (DockQ, CDR-H3 RMSD, epitope F1). The LLM translated this knowledge into steering objectives; without the knowledge, there would have been nothing to translate.

2. **Scientific literature**: The three publications provided the mathematical and conceptual scaffolding for novel approaches. The LLM's ability to read a preprint and synthesize an implementation — filling in unstated details, resolving ambiguities — was the critical bridge between published science and working code.

3. **An agentic coding LLM**: Not a chatbot answering questions, but a system that could read codebases, write and edit code, run evaluations, and interpret numerical results — completing the full engineering loop with minimal hand-holding.

---

## 3. The Broader Opportunity: Drug Discovery at Large

The antibody–antigen co-folding proof of concept is one instantiation of a pattern that applies across the field of computational drug discovery. The same three-ingredient recipe — domain expertise, recent literature, and an LLM coding agent — can be applied wherever deep learning models are used.

### 3.1 Candidate Application Areas

**3D Small Molecule Generative Models**  
Models such as DiffSBDD, Pocket2Mol, TargetDiff, and their successors generate drug-like molecules conditioned on protein binding sites. They are state-of-the-art for structure-based drug design (SBDD), but their outputs are often not optimized for the specific constraints of a real program — desired pharmacophores, obligate pharmacological features, ADMET liabilities. Inference-time steering techniques (exact analogues to what was done here) could guide these models toward regions of chemical space that satisfy program-specific multi-parameter optimization (MPO) criteria, without retraining.

**ADMET Property Prediction and Design**  
Transformer-based models for predicting absorption, distribution, metabolism, excretion, and toxicity (ADMET) properties are standard tools, but adapting them to program-specific liabilities — unusual scaffolds, rare metabolic pathways, species-specific toxicology data — typically requires either fine-tuning (expensive, data-hungry) or manual post-hoc filtering (wasteful). An LLM agent can implement soft-constraint steering techniques that incorporate ADMET objectives directly into the generative loop, using the same latent-space modulation approach demonstrated here.

**Protein Language Models for Sequence Design**  
ESM-2, ESMFold, and their successors are increasingly used for antibody humanization, CDR optimization, and affinity maturation. Inference-time guidance for these models — directing sequence generation toward desired biophysical properties, reduced immunogenicity flags, or specific developability profiles — is an active research area with several preprints published but few implementations. An LLM agent can translate these preprints into working implementations in the time it currently takes to review them.

**Molecular Dynamics and Free Energy Calculations**  
Neural network potentials (ANI, MACE, NequIP) and ML-enhanced free energy perturbation pipelines are increasingly used for binding affinity estimation. Customizing these pipelines to program-specific requirements — particular force field corrections, specialized sampling protocols for flexible binding sites, non-equilibrium work estimators — is exactly the kind of complex, paper-described modification that LLM-assisted implementation accelerates most dramatically.

**Multi-Modal Integration (Structure + Sequence + Experimental)**  
The most ambitious frontier is integrating multiple data modalities — AlphaFold-predicted structures, cryo-EM density maps, HDX-MS epitope data, SPR binding kinetics, and sequence databases — into unified inference pipelines. The EmbedOpt paper demonstrated that embedding-space optimization can naturally accommodate multiple experimental constraint types simultaneously. An LLM-assisted program could build institution-specific multi-modal steering pipelines that are simply not available as off-the-shelf tools.

### 3.2 The Common Thread

What unifies all these applications is the same pattern identified in the proof of concept:

> **Published methods describe what to do. LLMs can do it.**

The bottleneck was never the idea — scientific publications provide a steady stream of ideas. The bottleneck was the engineering capacity to translate ideas into production experiments. That bottleneck is now dramatically reduced.

---

## 4. Strategic Value and Competitive Advantage

### 4.1 Speed as a Competitive Moat

Drug discovery is a race. Decisions about which targets to pursue, which scaffolds to optimize, and which candidates to advance are made under uncertainty. The institutions that can most rapidly test a new method, generate a new hypothesis, or validate a new approach accumulate compounding advantages.

In the proof of concept, **27 significant improvements to existing strategies were designed and documented in a single analysis session** — based on reading three preprints and cross-referencing with the existing codebase. Two rounds of experiments were completed and fully analyzed. The cycle that typically spans a graduate student's thesis was compressed into weeks.

This speed does not come at the cost of rigor. Every strategy was benchmarked against ground-truth crystal structures with quantitative metrics and statistical significance testing. The failures were documented as carefully as the successes. The output is not a promising idea — it is a body of evaluated evidence.

### 4.2 Accessibility Across Departments

Historically, computational deep learning in drug discovery has been siloed in specialized groups. Medicinal chemists, structural biologists, PK/PD scientists, and toxicologists have not been able to directly customize ML models for their programs — they could only request help from the computational group, wait in the queue, and receive results that may or may not address their actual question.

LLM-assisted development changes this. A scientist who understands their problem well enough to describe it — the relevant biology, the key metrics, the constraints that matter — can now be a direct participant in model development. The LLM becomes an equalizer: it does not replace domain expertise, but it no longer requires a separate ML engineering expertise to act on domain knowledge.

This has organizational implications. Programs could run their own computational experiments, explore their own ideas, and generate their own evidence — with the speed and quality that previously required dedicated ML engineering resources.

### 4.3 Knowledge Accumulation

Each LLM-assisted experimental cycle produces structured, searchable documentation: implementation guides, evaluation reports, strategy comparisons, lessons-learned analyses. This is in contrast to traditional ML engineering work, where implementation knowledge is often tacit, undocumented, or lost when personnel change.

The proof of concept generated a corpus of documents — steering method analyses, round-by-round evaluations, improvement proposals — that constitutes institutional memory about what has been tried, what worked, and why. Future programs can build on this foundation rather than rediscovering the same lessons.

---

## 5. The Risk of Inaction

### 5.1 Competitive Erosion

The techniques described here are not proprietary. The papers are public. The LLMs are commercially available. Other institutions are already using them. The question is not whether this capability will exist in the field — it is whether our institution will be using it or observing it from behind.

The drug discovery industry is undergoing rapid AI adoption. Companies that establish LLM-assisted development workflows now will have:
- More experiments per unit time
- More diverse approaches tested
- More institutional knowledge accumulated
- Lower barriers to exploring adjacent scientific territory

Companies that wait will face a compounding disadvantage.

### 5.2 Talent and Resource Misallocation

The current model — bottlenecked by ML engineering capacity — systematically underutilizes domain expertise. Scientists who know what question to ask cannot pursue it without waiting for an engineer to translate it. Engineers spend time on implementation mechanics that an LLM can handle, rather than on the architectural decisions and experimental design that require genuine expertise.

Continuing this model means accepting that most of the domain knowledge in our institution will remain untested — not because the ideas aren't worth testing, but because the bandwidth to test them doesn't exist. LLM-assisted development is not about replacing scientists; it is about removing the engineering bottleneck that prevents scientists from doing science.

### 5.3 Publication and Priority

The proof of concept generated multiple strategies that represent genuine scientific contributions. The Hierarchical Steering approach (Strategy Y) — a two-phase inference method combining embedding-space guidance with coordinate-space refinement — achieved results that are, to our knowledge, unreported in the literature. Under a traditional workflow, the gap between "idea in a meeting" and "result worth publishing" would have been 12–18 months. Under LLM-assisted development, that gap was weeks.

Publication priority in computational drug discovery is increasingly determined by the speed of implementation cycles. Institutions that can move from paper → implementation → evaluation → publication faster than competitors will capture scientific credit for methods that others are thinking about but have not yet built.

---

## 6. Proposed Initiative

### 6.1 Phase 1: Expand the Proof of Concept (Months 1–3)

Extend the antibody–antigen steering work to address the key open problems identified in the two-round evaluation:

- Develop a custom re-ranking module to solve the confidence calibration bottleneck (the single highest-ROI intervention identified)
- Implement the Q→B2 pipeline: use LLM-assisted blind epitope discovery to generate contact restraints, removing the need for oracle structural information
- Test the combination of embedding-space steering (L+, Y) with coordinate-space contact restraints — a mechanistically orthogonal combination identified as high-promise

These are concrete, scoped experiments with clear success criteria, building directly on the existing benchmark.

### 6.2 Phase 2: Extend to a Second Modality (Months 3–6)

Apply the same LLM-assisted development workflow to a second model class — for example, a 3D small molecule generative model (DiffSBDD or TargetDiff). The objective is to demonstrate that the workflow is not specific to co-folding models but is a generalizable methodology.

Deliverables:
- Program-specific steering strategies for at least one active discovery program
- Benchmarked against relevant chemical property objectives (not just structural accuracy)
- Documentation suitable for handoff to the broader computational chemistry team

### 6.3 Phase 3: Institutionalize the Workflow (Months 6–12)

Establish the LLM-assisted development cycle as a standard operational mode for computational drug discovery at the institution:

- Develop reusable templates for strategy proposal, implementation, and evaluation
- Train scientists across departments to provide domain expertise as LLM inputs
- Establish a shared benchmark library for model evaluation
- Create a review process for incorporating insights from new publications into active programs

### 6.4 Success Metrics

| Phase | Metric | Target |
|-------|--------|--------|
| Phase 1 | DockQ improvement over vanilla baseline | > 0.05 without oracle information |
| Phase 1 | Time from new paper to evaluated implementation | < 2 weeks |
| Phase 2 | MPO score improvement in active program | Program-defined threshold |
| Phase 3 | Number of departments running own computational experiments | ≥ 3 |
| Phase 3 | Strategies evaluated per year | ≥ 20 |

---

## 7. Resource Requirements

The proof of concept was conducted with minimal dedicated resources: one scientist providing domain guidance, Claude Code as the LLM agent, and computational infrastructure already available for Boltz-2 predictions (GPU cluster). The primary resource constraint was not compute or personnel — it was the time required to design and write structured evaluation documents that could serve as LLM inputs.

For a full institutional program, the incremental requirements are:

- **Scientific coordinator** (0.5 FTE): Responsible for translating program questions into structured LLM inputs, reviewing outputs, and managing the feedback loop
- **Computational infrastructure**: Existing GPU resources are sufficient for Phase 1; Phase 2 and 3 may require modest additional capacity depending on the scale of benchmarking
- **LLM access**: Enterprise API or subscription access for agentic coding use
- **Literature monitoring**: A lightweight process for flagging relevant preprints to active programs — the raw material that feeds the implementation cycle

The expected return on these resources is measured in experiments that would otherwise not have been run, approaches that would otherwise have remained as ideas, and publications that would otherwise not have been written.

---

## 8. Conclusion

The proof of concept demonstrated that the barrier between scientific idea and working implementation has dropped dramatically. A co-folding model was extended with novel inference-time steering strategies derived directly from preprints — strategies that would previously have required months of specialized ML engineering work to implement and evaluate. Two complete experimental rounds were completed, 17 strategies were benchmarked, and the results advanced the state of the art in blind antibody–antigen docking without oracle information.

The generalization of this capability to the full breadth of drug discovery tools — generative models for small molecules, sequence-based antibody design, ML-enhanced free energy calculations, multi-modal data integration — is not speculative. It follows directly from the same methodology, applied to different scientific domains.

The central message is this: **improving and customizing deep learning models for drug discovery is more accessible than ever before — to every department, every program, and every scientist who can articulate what a good result looks like.** The institutions that embrace this new paradigm will move faster, explore more ideas, and accumulate more knowledge than those that do not. The window for establishing an early advantage is open now, but it will not remain open indefinitely.

We recommend approving this initiative and beginning Phase 1 immediately.

---

## Appendix: Proof-of-Concept Summary Statistics

| Category | Value |
|----------|-------|
| Model | Boltz-2 (state-of-the-art antibody–antigen co-folding model) |
| Benchmark | 47 antibody–antigen complexes, ground-truth crystal structures |
| Experimental rounds | 2 (Round 1: 11 strategies; Round 2: 6 strategies) |
| Source papers incorporated | 3 (FK-Diffusion Steering; Boltz-sample; EmbedOpt) |
| Novel strategies proposed and implemented | 17 |
| Improvements and enhancements documented | 27 |
| Best DockQ improvement (Strategy A) | +23% over baseline (0.247 vs 0.201) |
| Statistically significant result | CDR-H3 loop error –0.30 Å (Strategy Y, p = 0.025) |
| Best blind epitope prediction improvement | +60% relative gain (Strategy Q, median F1) |
| Time from paper to evaluated implementation | Days to weeks (vs. typical months) |
| Oracle vs. confidence selection gap (A+) | 0.044 DockQ — generation quality already exceeds oracle target; selection is the remaining bottleneck |
