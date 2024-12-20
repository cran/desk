## desk 1.1.2

# Additions:
-   package vignette

# Fixes:
-   year of publication of `Ökonometrie - Das R-Arbeitsbuch, 2nd ed.`: 2023 to 2024
-   fix in `ols` when decimal separator is switched to a comma

## desk 1.1.2

# Fixes:
-   Suggested packages imported conditionally in examples

## desk 1.1.1

# Modifications:
-   Change in `plot.desk`: legend without background color
-   Minor changes in `data.eu`, `bp.test`, `roll.win`, `wh.test`
-   Change variable name `exp` to `expend` in `data.govexpend`
-   Change variable names `f_wage`, `a_wage`, `f_age`, `a_age` to 
    `f.wage`, `a.wage`, `f.age`, `a.age` in `data.windscreen`
    
# Additions:
-   New variable `religion` in data.wage

# Fixes:
-   Minor fix in `new.session`

## desk 1.1.0

This is the first official CRAN release

# Modifications (compared to non-CRAN 1st ed. releases):
- Function renames
    `def.exp` (old: exp.def)
    `def.log` (old: log.def)
    `par.f.test` (old: f.coef.test)
    `par.t.test` (old: t.coef.test)
    `repeat.sample` (old: rep.sample)

# Additions:
-   New functions
    `new.session()` 
    `roll.win()`
-   New datasets
    `data.macro` 
    `data.regional` 
    `data.spurious`
-   New welcome screen with download link for 1.0.x version (first edition of book)

# Fixes:
-   dataset and function help files updated and completed
