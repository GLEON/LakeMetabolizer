## Test environments
- local macOS, R 4.5.0
- macOS (GitHub Actions, `macos-latest`), R release
- Windows (GitHub Actions, `windows-latest`), R release
- Ubuntu (GitHub Actions, `ubuntu-latest`), R devel (with release UA), _release_, _oldrel-1_

## Continuous integration details
- R CMD check via `r-lib/actions/check-r-package@v2`, which runs `rcmdcheck` with `--as-cran` on each matrix entry.
- Dependency setup with `r-lib/actions/setup-r-dependencies@v2`; caches deps and installs `rcmdcheck` and `covr`.
- Pandoc available (`setup-pandoc@v2`) so vignettes and Rmd docs can build in CI.
- Code coverage computed with `covr::codecov()` after checks on each job; reports sent to Codecov.

## R CMD check results
0 errors | 0 warnings | 1 notes
- NOTE: standardizing maintainer name to `Jacob Zwart`; previously `Jake Zwart`

* This is a new release. 1.5.6 

## Notes for CRAN
Maintenance release: modernize inst/CITATION, fix Rd cross-refs, standardize maintainer name, regenerate docs. No functional changes.

