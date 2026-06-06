# go-kriging

Flat Go library (`package kriging`) for Kriging spatial interpolation with GeoTIFF raster output.

## Build / Test

Standard Go tooling — no Makefile, no CI, no pre-commit hooks.

```bash
go build ./...
go test ./...
go vet ./...
go fmt ./...
```

Run a single test:
```bash
go test -v -run TestNewConvex ./...
go test -v -run TestInterpolator ./...
```

## Local dependencies (gotcha)

`go.mod` has two `replace` directives pointing to sibling directories:

```
replace github.com/flywave/go-geos => ../go-geos
replace github.com/flywave/go-geoid => ../go-geoid
```

**`go build` / `go test` will fail** unless `../go-geos` and `../go-geoid` are present. These are not published remotely — they must exist locally.

## Package structure

All source is in the root package `kriging` — no sub-packages.

| File | Role |
|------|------|
| `kriging.go` | Core: `Train()`, `Predict()`, `Contour()` |
| `interpolator.go` | Composite: `KrigingInterpolator`, `BilinearInterpolator`, `HyperbolicInterpolator` |
| `convex.go` | QuickHull 2D convex hull |
| `voxelgrid.go` | Point cloud downsampling via voxel grid filter |
| `grid.go` | Raster grid generation |
| `matrix.go` / `matrix-inverse.go` | Matrix ops (transpose, multiply, Cholesky, solve, inverse via gonum fallback) |

## Test data (committed)

Integration tests read real files from the repo root:

- `test.json` — GeoJSON FeatureCollection (manhole elevation points)
- `02.geojson` — GeoJSON FeatureCollection (building footprint points)
- `test.tif` — background DEM GeoTIFF

Output GeoTIFFs (`out1.tif`, `out2.tif`) are also committed — treat as generated artifacts.

## Variogram models

`gaussian`, `exponential`, `spherical` — passed as `kriging.ModelType` strings.

## Framework

- Tests use `github.com/stretchr/testify/assert` — same package, not `_test` package.
- Uses `ioutil.ReadAll` (deprecated since Go 1.16); prefer `io.ReadAll` in new code.
