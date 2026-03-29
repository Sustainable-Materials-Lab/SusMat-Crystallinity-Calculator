## [2.0.5](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v2.0.4...v2.0.5) (2026-03-29)


### Bug Fixes

* add pypi release via Github ([ebb4673](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/ebb4673b059265419e2c3b5313f2201dffd35d45))

## [2.0.4](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v2.0.3...v2.0.4) (2026-02-24)


### Bug Fixes

* add dependencies to avoid vulnerable matplotlib dependencies. ([d759ca5](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/d759ca5a79a6a2e8ef3b64875b5e274b50581899))

## [2.0.3](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v2.0.2...v2.0.3) (2026-02-24)


### Bug Fixes

* add debug print to data importer ([62ad622](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/62ad6226f80fa7dd18824a7d25c2dd94f659366f))

## [2.0.2](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v2.0.1...v2.0.2) (2026-02-13)


### Bug Fixes

* typo preventing recognition of error column. Fixes [#3](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/issues/3) ([ddc5d2d](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/ddc5d2da9cff1fa474cc61435361469d78de94c8))

## [2.0.1](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v2.0.0...v2.0.1) (2025-11-12)


### Bug Fixes

* enable vuln check with trivy ([93aea72](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/93aea7263e2ccd990f5a3df60c00808d346eb24d))

# [2.0.0](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v1.0.0...v2.0.0) (2025-11-03)


### Features

* **core:** switch to UV build system ([3db31dd](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/3db31dd8c73f7302137a9dc29af21bf9003030f9))


### BREAKING CHANGES

* **core:** switch build system

# [1.0.0](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v0.3.0...v1.0.0) (2025-10-31)


### Features

* added packaging metadata ([5e8d10d](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/5e8d10dfebf9b29487fd31d077afc2879a12382c))


### BREAKING CHANGES

* no longer namespace package

# [0.3.0](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v0.2.1...v0.3.0) (2025-10-30)


### Features

* **core:** autodetect presence of standard errors in input data ([9f54c1b](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/9f54c1b5b888a1fe04582c8e6923e4cd9f9a46d1))

## [0.2.1](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v0.2.0...v0.2.1) (2025-06-16)


### Bug Fixes

* **remote:** explicitly state output format for saving figures rather than relying on file extension. Fixes [#1](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/issues/1) ([a96226b](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/a96226b00df3dca68f93b5133bcef81578675175))

# [0.2.0](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v0.1.2...v0.2.0) (2025-03-09)


### Bug Fixes

* legend not showing Jeffamine label ([d51b193](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/d51b193408c11e29610fff7083c6ce2d8cfca63d))


### Features

* detect and import jeffamine data for local files ([1e04353](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/1e04353255db049ab4450397870b88fca4ba2f5f))

## [0.1.2](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/compare/v0.1.1...v0.1.2) (2025-03-05)


### Bug Fixes

* match crystalline integral curve colour to the colour of the cellulose I curve ([d7164a7](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/d7164a7815ab243a2f2c957bf10b79161333501b))
* removed scattering vector unit conversion to nm (consistent with XSACT and other software output). ([8b20afe](https://gitlab.kuleuven.be/susmat/cellulose/waxs-crystallinity-calculator/commit/8b20afe55846616df60a4923a3e7ecd5d1ee0f09))
