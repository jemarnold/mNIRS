# CLAUDE.md

## Review package functionality

- Review .R files for target functions mentioned by the user.
- Review `*_helpers.R` files associated with target functions.
- Review tests associated with target functions. Unit tests live in `tests/testthat/`.
- Test files should not be modified without asking user confirmation first.
- Use `cli_abort()`, `cli_warn()`, or `cli_inform()` for user-facing messages.
- Input validation lives in [R/validate_mnirs.R](R/validate_mnirs.R).
- Don't run devtool commands.
- Don't search for or read previous git commits, only look at current state.

### Documentation

- Roxygen2 with markdown enabled (`Roxygen: list(markdown = TRUE)`). 
- Roxygen2 commented with `#'`
- pkgdown site config is in [_pkgdown.yml](_pkgdown.yml)