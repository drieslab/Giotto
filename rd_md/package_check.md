# `package_check`

package_check


## Description

check if package is available and provide installation instruction if not available


## Usage

```r
package_check(
  pkg_name,
  repository = c("CRAN", "Bioc", "github", "pip"),
  github_repo = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`pkg_name`     |     name of package
`repository`     |     where is the package
`github_repo`     |     name of github repository if needed


