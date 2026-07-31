## Review package functionality

- Review AGENTS files in `/inst/skills/` for package functionality overview.
- Unit tests live in `tests/testthat/`.
- Input validation lives in [R/validate_mnirs.R](R/validate_mnirs.R).
- Roxygen2 with markdown enabled (`Roxygen: list(markdown = TRUE)`) commented with `#'`.
- pkgdown site config is in [_pkgdown.yml](_pkgdown.yml)

## Scope limitations
- Don't run tests, build, checks, render documentation, or any other devtool/package commands through `devtool` or `Rscript`.
- Don't make any write changes to git/GitHub, read only access.
- Don't search for or read previous git commits, only look at current state.
- Don't manually edit Roxygen2-generated .Rd files.
- Package documentation (NEWS.md, DESCRIPTION, NAMESPACE, README, vignettes, articles etc.) should not be modified without asking user confirmation first.
- Test files should not be modified without asking user confirmation first.