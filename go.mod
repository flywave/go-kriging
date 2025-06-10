module github.com/flywave/go-kriging

go 1.23.0

toolchain go1.24.2

require (
	github.com/flywave/go-cog v0.0.0-20250607133043-41acd04eb904
	github.com/flywave/go-geo v0.0.0-20250607132733-46bd30e585ce
	github.com/flywave/go-geoid v0.0.0-20210705014121-cd8f70cb88bb
	github.com/flywave/go-geom v0.0.0-20250607125323-f685bf20f12c
	github.com/flywave/go3d v0.0.0-20250314015505-bf0fda02e242
	github.com/stretchr/testify v1.10.0
	gonum.org/v1/gonum v0.8.2
)

require (
	github.com/davecgh/go-spew v1.1.1 // indirect
	github.com/flywave/go-geos v0.0.0-20250607125930-047054a9f657 // indirect
	github.com/flywave/go-proj v0.0.0-20250607132305-d70d32f5ad2d // indirect
	github.com/google/tiff v0.0.0-20161109161721-4b31f3041d9a // indirect
	github.com/hhrutter/lzw v1.0.0 // indirect
	github.com/kr/text v0.2.0 // indirect
	github.com/pmezard/go-difflib v1.0.0 // indirect
	golang.org/x/exp v0.0.0-20191030013958-a1ab85dbe136 // indirect
	golang.org/x/image v0.28.0 // indirect
	gopkg.in/yaml.v3 v3.0.1 // indirect
)

replace github.com/flywave/go-geos => ../go-geos

replace github.com/flywave/go-geoid => ../go-geoid
