package kriging

import (
	"math"
	"testing"

	vec3d "github.com/flywave/go3d/float64/vec3"
	"github.com/stretchr/testify/assert"
)

func TestMinMaxFloat64(t *testing.T) {
	assert := assert.New(t)

	allPos := []vec3d.T{{10, 20, 0}, {30, 50, 0}, {5, 100, 0}}
	assert.Equal(float64(5), minFloat64(allPos, 0))
	assert.Equal(float64(20), minFloat64(allPos, 1))
	assert.Equal(float64(30), maxFloat64(allPos, 0))
	assert.Equal(float64(100), maxFloat64(allPos, 1))

	allNeg := []vec3d.T{{-10, -20, 0}, {-30, -5, 0}, {-1, -100, 0}}
	assert.Equal(float64(-30), minFloat64(allNeg, 0))
	assert.Equal(float64(-100), minFloat64(allNeg, 1))
	assert.Equal(float64(-1), maxFloat64(allNeg, 0))
	assert.Equal(float64(-5), maxFloat64(allNeg, 1))

	mixed := []vec3d.T{{-10, 20, 0}, {30, -50, 0}, {0, 0, 0}}
	assert.Equal(float64(-10), minFloat64(mixed, 0))
	assert.Equal(float64(-50), minFloat64(mixed, 1))
	assert.Equal(float64(30), maxFloat64(mixed, 0))
	assert.Equal(float64(20), maxFloat64(mixed, 1))

	empty := []vec3d.T{}
	assert.Equal(float64(0), minFloat64(empty, 0))
	assert.Equal(float64(0), maxFloat64(empty, 0))
}

func TestContourAxisSwap(t *testing.T) {
	assert := assert.New(t)

	points := []vec3d.T{
		{0, 0, 10},
		{20, 0, 20},
		{0, 10, 30},
		{20, 10, 40},
	}
	kri := New(points)
	_, err := kri.Train(Spherical, 0, 100)
	if err != nil {
		t.Fatal(err)
	}

	xWidth, yWidth := 2, 3
	rect := kri.Contour(xWidth, yWidth)

	assert.Equal(xWidth, rect.XWidth)
	assert.Equal(yWidth, rect.YWidth)
	assert.Equal(float64(0), rect.Xlim[0])
	assert.Equal(float64(20), rect.Xlim[1])
	assert.Equal(float64(0), rect.Ylim[0])
	assert.Equal(float64(10), rect.Ylim[1])

	xl := rect.Xlim[1] - rect.Xlim[0]
	yl := rect.Ylim[1] - rect.Ylim[0]
	gridW := xl / float64(xWidth)
	gridH := yl / float64(yWidth)

	for j := 0; j < yWidth; j++ {
		yTarget := rect.Ylim[0] + float64(j)*gridH
		for k := 0; k < xWidth; k++ {
			xTarget := rect.Xlim[0] + float64(k)*gridW
			idx := j*xWidth + k
			expected := kri.Predict(xTarget, yTarget)
			assert.InDelta(expected, rect.Contour[idx], 1e-6,
				"contour[%d] = %.6f, expected Predict(%.2f, %.2f) = %.6f",
				idx, rect.Contour[idx], xTarget, yTarget, expected)
		}
	}
}

func TestVoxelFilterNoSkip(t *testing.T) {
	assert := assert.New(t)

	points := make([]vec3d.T, 10)
	for i := 0; i < 10; i++ {
		points[i] = vec3d.T{float64(i * 100), float64(i * 100), 0}
	}

	leafSize := vec3d.T{50, 50, 50}
	vg := newVoxelGrid(leafSize)

	result, err := vg.Filter(points)
	assert.NoError(err)

	assert.Len(result, 10)

	for i := 0; i < 10; i++ {
		assert.Equal(float64(i*100), result[i][0])
		assert.Equal(float64(i*100), result[i][1])
	}
}

func TestVoxelFilterSinglePointInVoxel(t *testing.T) {
	assert := assert.New(t)

	points := make([]vec3d.T, 4)
	points[0] = vec3d.T{1, 2, 0}
	points[1] = vec3d.T{2, 2, 0}
	points[2] = vec3d.T{100, 200, 0}
	points[3] = vec3d.T{150, 220, 0}

	leafSize := vec3d.T{50, 50, 50}
	vg := newVoxelGrid(leafSize)

	result, err := vg.Filter(points)
	assert.NoError(err)

	assert.Len(result, 3)

	foundCentroid := false
	foundFirstSingle := false
	foundSecondSingle := false
	for _, p := range result {
		if math.Abs(p[0]-1.5) < 1e-10 {
			foundCentroid = true
		}
		if math.Abs(p[0]-100) < 1e-10 {
			foundFirstSingle = true
		}
		if math.Abs(p[0]-150) < 1e-10 {
			foundSecondSingle = true
		}
	}
	assert.True(foundCentroid, "should have centroid of merged voxel at x=1.5")
	assert.True(foundFirstSingle, "single point at x=100 should retain its coordinate")
	assert.True(foundSecondSingle, "single point at x=150 should retain its coordinate")
}

func TestMatrixMultiplyZerosPreserved(t *testing.T) {
	assert := assert.New(t)

	a := []float64{
		0, 5,
		3, 0,
		0, 0,
	}
	b := []float64{
		1, 2,
		3, 4,
	}

	result := matrixMultiply(a, b, 3, 2, 2)

	expected := []float64{
		15, 20,
		3, 6,
		0, 0,
	}
	for i := range result {
		assert.Equal(expected[i], result[i])
	}
}

func TestMatrixCholAccuracy(t *testing.T) {
	assert := assert.New(t)

	a := []float64{
		4, 2, 1, 1,
		2, 5, 3, 2,
		1, 3, 6, 4,
		1, 2, 4, 7,
	}
	original := make([]float64, 16)
	copy(original, a)

	ok := matrixChol(a, 4)
	assert.True(ok)

	reconstructed := make([]float64, 16)
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			sum := 0.0
			for k := 0; k < 4; k++ {
				li_k := 0.0
				if k < i {
					li_k = a[i*4+k]
				} else if k == i {
					li_k = a[i*4+i]
				}
				lj_k := 0.0
				if k < j {
					lj_k = a[j*4+k]
				} else if k == j {
					lj_k = a[j*4+j]
				}
				sum += li_k * lj_k
			}
			reconstructed[i*4+j] = sum
		}
	}

	for i := range original {
		assert.InDelta(original[i], reconstructed[i], 1e-10,
			"A[%d] = %.10f, reconstructed = %.10f", i, original[i], reconstructed[i])
	}
}
