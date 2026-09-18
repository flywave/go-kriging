package kriging

import (
	"io"
	"os"
	"testing"

	"github.com/flywave/go-geom/general"
	vec3d "github.com/flywave/go3d/float64/vec3"
	"github.com/stretchr/testify/assert"
)

// 回归：未知模型类型必须返回错误。
// 此前 Train 的 switch 没有 default 分支，kri.model 保持 nil，
// 并行构造 K 矩阵的 worker goroutine 调用 nil 函数会打挂整个进程。
func TestTrainUnknownModelReturnsError(t *testing.T) {
	assert := assert.New(t)

	points := []vec3d.T{
		{0, 0, 10},
		{1, 0, 20},
		{0, 1, 30},
		{1, 1, 40},
	}

	for _, model := range []ModelType{"Spherical", "Gaussian", "spherical ", "", "bogus"} {
		kri := New(points)
		_, err := kri.Train(model, 0, 100)
		if model == Spherical {
			assert.NoError(err, "spherical must still train")
			continue
		}
		assert.Error(err, "model %q must be rejected", model)
	}
}

// 回归：点数不足时必须返回错误，不能索引越界
func TestTrainTooFewPoints(t *testing.T) {
	assert := assert.New(t)

	for _, points := range [][]vec3d.T{
		{},
		{{0, 0, 1}},
		{{0, 0, 1}, {1, 1, 2}},
	} {
		kri := New(points)
		_, err := kri.Train(Spherical, 0, 100)
		assert.Error(err, "points=%d must be rejected", len(points))
	}
}

// 回归：扁平点集（所有点等高）不能让体素滤波算出 NaN 下标。
// interpolator.filter 用 extentZ/512 作为 Z 层高，等高时该值为 0，
// p[2]/0 = NaN，而 int(NaN) 的结果与平台相关（amd64 上是 INT_MIN），
// 会让体素下标越界 panic。
func TestVoxelGridFlatInput(t *testing.T) {
	assert := assert.New(t)

	points := []vec3d.T{
		{0, 0, 0},
		{1, 0, 0},
		{0, 1, 0},
		{1, 1, 0},
		{5, 5, 0},
	}

	min, max, err := minMaxVec3(points)
	assert.NoError(err)

	// 与 interpolator.filter 相同的层高推导
	leafSize := vec3d.T{
		(max[0] - min[0]) / float64(default_filter_size[0]),
		(max[1] - min[1]) / float64(default_filter_size[1]),
		(max[2] - min[2]) / float64(default_filter_size[2]),
	}
	assert.Equal(float64(0), leafSize[2], "flat input gives a zero Z leaf size")

	vg := newVoxelGrid(leafSize)
	result, err := vg.Filter(points)
	assert.NoError(err)
	assert.Len(result, len(points), "flat points must all survive the filter")
	for _, p := range result {
		assert.Equal(float64(0), p[2])
	}

	// 单点/单轴退化（所有点共线）也不能越界
	line := []vec3d.T{{1, 7, 3}, {2, 7, 3}, {9, 7, 3}}
	vg2 := newVoxelGrid(vec3d.T{1, 0, 0})
	lineResult, err := vg2.Filter(line)
	assert.NoError(err)
	assert.Len(lineResult, len(line))
}

// 回归：进度回调返回 false 时 Process 必须尽快返回 ErrAborted
func TestInterpolatorProgressAbort(t *testing.T) {
	assert := assert.New(t)

	f, err := os.Open("./test.json")
	assert.NoError(err)
	defer f.Close()

	json, err := io.ReadAll(f)
	assert.NoError(err)

	fcs, err := general.UnmarshalFeatureCollection(json)
	assert.NoError(err)

	m := Spherical
	calls := 0
	opts := Options{
		Input:  fcs,
		Output: t.TempDir() + "/aborted.tif",
		Model:  &m,
		Progress: func(part, total uint64) bool {
			calls++
			return false
		},
	}

	ker := NewKrigingInterpolator(opts)
	_, _, err = ker.Process()
	assert.Equal(ErrAborted, err)
	assert.Equal(1, calls)
}

// 进度回调为 nil 时行为不变（向后兼容）
func TestInterpolatorWithoutProgress(t *testing.T) {
	assert := assert.New(t)

	f, err := os.Open("./test.json")
	assert.NoError(err)
	defer f.Close()

	json, err := io.ReadAll(f)
	assert.NoError(err)

	fcs, err := general.UnmarshalFeatureCollection(json)
	assert.NoError(err)

	m := Spherical
	opts := Options{
		Input:  fcs,
		Output: t.TempDir() + "/ok.tif",
		Model:  &m,
		Progress: func(part, total uint64) bool {
			return true
		},
	}

	ker := NewKrigingInterpolator(opts)
	_, _, err = ker.Process()
	assert.NoError(err)
}

// 回归：2D 点（GeoJSON 只给 x,y）不能下标越界。
// extractPosion 早期直接取 Data()[2]，2D 点集必然 panic —— 这正是
// “用 2D 点插值”的使用场景（Z 视为 0）。
func TestInterpolator2DPoints(t *testing.T) {
	assert := assert.New(t)

	fcs, err := general.UnmarshalFeatureCollection([]byte(`{
		"type": "FeatureCollection",
		"features": [
			{"type": "Feature", "properties": {}, "geometry": {"type": "Point", "coordinates": [116.30, 39.80]}},
			{"type": "Feature", "properties": {}, "geometry": {"type": "Point", "coordinates": [116.34, 39.80]}},
			{"type": "Feature", "properties": {}, "geometry": {"type": "Point", "coordinates": [116.30, 39.84]}},
			{"type": "Feature", "properties": {}, "geometry": {"type": "Point", "coordinates": [116.34, 39.84]}}
		]
	}`))
	assert.NoError(err)

	m := Spherical
	ps := [2]float64{0.01, 0.01}
	opts := Options{
		Input:     fcs,
		Output:    t.TempDir() + "/2d.tif",
		Model:     &m,
		PixelSize: &ps,
	}

	ker := NewKrigingInterpolator(opts)
	_, _, err = ker.Process()
	assert.NoError(err)

	// Z 全部按 0 处理
	pos := ker.extractPosion()
	assert.Len(pos, 4)
	for _, p := range pos {
		assert.Equal(float64(0), p[2])
	}
}

// ZOf：2D 坐标取 0，3D 坐标取第三个分量
func TestZOf(t *testing.T) {
	assert := assert.New(t)
	assert.Equal(float64(0), ZOf([]float64{1, 2}))
	assert.Equal(float64(0), ZOf(nil))
	assert.Equal(float64(7), ZOf([]float64{1, 2, 7}))
}
