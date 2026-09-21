# Final submission-readiness audit

## Status

**Scientific and numerical consistency: PASS after final corrections.**

The compact manuscript source package was audited after Tasks 1-8 for numerical consistency, stale claims, citation/reference integrity, build reproducibility, package documentation and PDF layout. A repeatable compact-source check is included as `tests/test_submission_consistency.py`.

## Final corrections

- Regenerated the BibTeX bibliography so numeric citations follow first-citation order after the Task 8 related-work additions.
- Tightened the description of scDisInFact to distinguish shared/unshared biological factors from explicit batch covariates and alignment.
- Corrected the data/code availability statement: the compact manuscript package contains manuscript tables, figures, model summaries, gene lists and figure-generation code, but not the large dense arrays, exact per-cell split-index archive or fitted checkpoint binaries.
- Replaced the compact-package root README with bundle-specific build and reproduction instructions.
- Updated stale audit documentation, including the canonical cell-level split versus the separate Task 5 complete-library holdout and removal of obsolete author-asterisk language.
- Corrected the paper-table generator so regeneration preserves the scDisInFact software-version row and no longer emits superseded historical pathway macros.
- Added PDF title/author metadata and rebuilt the manuscript and supplement from source.

## Verification

The final compact-source audit checks cohort totals, split sizes, Task 2-7 headline results, citation keys, LaTeX labels, stale permutation claims, pathway-macro removal and package-availability wording. The final build has no undefined citations, undefined references, overfull boxes or fatal LaTeX errors.

Final page counts:
- main manuscript: 13 pages;
- supplementary material: 12 pages.

## Remaining submission-specific items

The scientific package is internally consistent, but journal submission is not administratively complete until the authors supply or confirm:

- corresponding author and contact information;
- author-contribution statement;
- funding statement;
- competing-interest statement;
- final data-reuse/data-availability wording;
- target journal, article type and journal-specific formatting/reporting checklist;
- any journal-specific animal-ethics wording appropriate for computational reanalysis of public animal data.

If definitive animal-level sample totals are reported later, the source-metadata pooling details for GSE216146 should be reconciled explicitly rather than inferred from deposited library counts.

An immutable repository release or DOI is also preferable before submission.
