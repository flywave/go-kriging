# go-kriging

Go library for Kriging spatial interpolation with GeoTIFF raster output.

## Features

- **Kriging interpolation** — Gaussian, Exponential, Spherical variogram models
- **Composite interpolation** — uses Kriging inside the convex hull, Bilinear/Hyperbolic resampling of a background DEM outside
- **GeoTIFF output** — writes Cloud-Optimized GeoTIFF (COG) via `go-cog`
- **Voxel grid filter** — point cloud downsampling before Kriging training
- **Convex hull** — QuickHull 2D algorithm, used to delimit the interpolation region

## Usage

```go
import "github.com/flywave/go-kriging"
```

### Core Kriging

```go
points := []vec3d.T{{lon, lat, elev}, ...}
krig := kriging.New(points)
krig.Train(kriging.Spherical, 0, 100)

z := krig.Predict(lon, lat)          // predict at a point
rect := krig.Contour(width, height)   // generate a prediction grid
```

### Composite interpolator (Kriging + background DEM)

```go
opts := kriging.Options{
    Input:    featureCollection,   // *geom.FeatureCollection
    Output:   "./out.tif",
    Model:    modelType,           // kriging.Spherical / Gaussian / Exponential
}
inter := kriging.NewKrigingInterpolator(opts)
inter.Process()                     // train, resample, write GeoTIFF
```

See `interpolator_test.go` for full integration examples.

## Commands

```bash
go build ./...
go test ./...
go vet ./...
go fmt ./...
```

### Single test

```bash
go test -v -run TestNewConvex ./...
go test -v -run TestInterpolator1 ./...
```

## Variogram models

| Constant     | String         |
|--------------|----------------|
| `Gaussian`   | `"gaussian"`   |
| `Exponential`| `"exponential"`|
| `Spherical`  | `"spherical"`  |

## Package structure

All source is in the root package `kriging` — no sub-packages.

| File | Role |
|------|------|
| `kriging.go` | `Train()`, `Predict()`, `Contour()` |
| `interpolator.go` | `KrigingInterpolator`, `BilinearInterpolator`, `HyperbolicInterpolator` |
| `convex.go` | QuickHull 2D convex hull |
| `voxelgrid.go` | Voxel grid point cloud filter |
| `grid.go` | Raster grid generation |
| `matrix.go` / `matrix-inverse.go` | Matrix ops (transpose, multiply, Cholesky, solve) |

## Dependencies (gotcha)

`go.mod` has two `replace` directives pointing to sibling directories:

```text
replace github.com/flywave/go-geos => ../go-geos
replace github.com/flywave/go-geoid => ../go-geoid
```

`go build` / `go test` will fail unless `../go-geos` and `../go-geoid` are present locally. These are not published remotely.

Key direct dependencies: `go-cog`, `go-geo`, `go-geoid`, `go-geom`, `go3d`, `gonum`.

## Test data

Files committed in the repo root:

- `test.json` / `02.geojson` — GeoJSON FeatureCollection inputs
- `test.tif` — background DEM GeoTIFF

Generated output `out1.tif` / `out2.tif` / `out21.tif` are committed as artifacts.

## Performance

`TestInterpolator1` (~1400 training points, full Kriging training + grid resampling + GeoTIFF write) runs in ~10s on Apple M4. The O(n³) Cholesky decomposition uses gonum's LAPACK backend.
