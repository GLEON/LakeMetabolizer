# Version 1.5.6 (2025-10-30)

-   Modernized `inst/CITATION` to current `bibentry()` style to resolve CRAN NOTE about old-style citation entries.
-   Fixed Rd cross-references by adding package-qualified anchors for functions in `rLakeAnalyzer`.
-   Regenerated documentation with roxygen2 7.3.3; updated `RoxygenNote` and rebuilt `NAMESPACE`.
-   Refreshed package-level help (`R/LakeMetabolizer-package.R`) with updated author/maintainer text.
-   DESCRIPTION maintenance: bumped version to 1.5.6 and standardized Maintainer name.
-   Other minor housekeeping

# Version 1.5.5 (2022-10-29)

-   Changes default branch on GitHub from `master` to `main`

-   Migrated away from travis-ci to GitHub Actions for package builds and checks

-   Updated recommended citation

-   Added codecov CI

# Version 1.5.3 (2020-04-23)

-   Fixes `k.read`, `k.read.soloviev`, and `k.macIntyre` bug where base functions fail when calculating kinematic viscosity

# Version 1.5.1 (2017-02-12)

-   Fix `k.heiskanen` bug where it fails when trying to calculate lwnet
