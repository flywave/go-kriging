package kriging

import (
	"io"
	"math"
	"os"
	"strings"
	"testing"

	"github.com/flywave/go-geom/general"
	vec3d "github.com/flywave/go3d/float64/vec3"
	"github.com/stretchr/testify/assert"
)

// 病态系统必须显式报错：gaussian 在 test.json 这种小 range 大点集上
// 条件数 ~1e20，旧代码会解出 1e45 量级的权重、输出 1e30 量级的预测值。
func TestTrainIllConditionedError(t *testing.T) {
	assert := assert.New(t)

	pos := loadVerifyPoints(t)
	kri, err := New(pos).Train(Gaussian, 0, 100)
	if err == nil {
		// 若恰好可解，结果也必须是有限的
		for i := range kri.pos {
			v := kri.Predict(kri.pos[i][0], kri.pos[i][1])
			assert.False(math.IsNaN(v) || math.IsInf(v, 0))
		}
		return
	}
	assert.Nil(kri)
	// gonum 直接判奇异，或我们的残差检查报病态，二者皆可
	if !strings.Contains(err.Error(), "singular") {
		assert.Contains(err.Error(), "ill-conditioned")
	}
}

// spherical 在同一数据上必须精确插值（sigma2=0）
func TestTrainExactAtDataPoints(t *testing.T) {
	assert := assert.New(t)

	pos := loadVerifyPoints(t)
	kri, err := New(pos).Train(Spherical, 0, 100)
	assert.NoError(err)
	for i := range kri.pos {
		e := math.Abs(kri.Predict(kri.pos[i][0], kri.pos[i][1]) - kri.pos[i][2])
		assert.Less(e, 1e-6, "point %d", i)
	}
}

// z 全等的平面点集：拟合斜率为 0，K 全零。必须降级为常数（均值）
// 预测器而不是奇异失败或静默输出垃圾。
func TestTrainFlatZMeanPredictor(t *testing.T) {
	assert := assert.New(t)

	pos := []vec3d.T{
		{116.30, 39.80, 7}, {116.34, 39.80, 7}, {116.30, 39.84, 7}, {116.34, 39.84, 7},
	}
	kri, err := New(pos).Train(Spherical, 0, 100)
	assert.NoError(err)
	for _, p := range [][2]float64{{116.31, 39.81}, {116.33, 39.83}} {
		assert.InDelta(7.0, kri.Predict(p[0], p[1]), 1e-12)
	}
}

// 拟合可能给出负块金（无约束回归），必须钳制为非负
func TestTrainNuggetClamped(t *testing.T) {
	pos := loadVerifyPoints(t)
	kri, err := New(pos).Train(Spherical, 0, 100)
	if err != nil {
		t.Skip("train failed on this dataset")
	}
	assert.GreaterOrEqual(t, kri.nugget, 0.0)
}

// 平面距离全部相同的点集（等边三角形）：range 退化接近 0，
// 不能产生 NaN/Inf
func TestTrainDegenerateRange(t *testing.T) {
	assert := assert.New(t)

	pos := []vec3d.T{
		{0, 0, 1}, {1, 0, 2}, {0.5, math.Sqrt(3) / 2, 3},
	}
	kri, err := New(pos).Train(Spherical, 0, 100)
	if err != nil {
		return // 显式报错也可接受
	}
	v := kri.Predict(0.3, 0.3)
	assert.False(math.IsNaN(v) || math.IsInf(v, 0))
}

// 近重合点必须被合并：否则 K 矩阵到求解器精度下奇异
func TestTrainDedupNearCoincident(t *testing.T) {
	assert := assert.New(t)

	pos := []vec3d.T{
		{0, 0, 10}, {0, 0, 12}, // 完全重合
		{1e-9, 1e-9, 11}, // 近重合
		{1, 0, 20},
		{0, 1, 30},
		{1, 1, 40},
	}
	kri, err := New(pos).Train(Spherical, 0, 100)
	assert.NoError(err)
	assert.Less(len(kri.pos), len(pos))
	// 合并后求解必须成功且预测有限（该玩具数据斜率可能为 0，
	// 命中均值预测器退化为常数也合法）
	for i := range kri.pos {
		v := kri.Predict(kri.pos[i][0], kri.pos[i][1])
		assert.False(math.IsNaN(v) || math.IsInf(v, 0), "point %d = %v", i, v)
	}
}

// 体素滤波不能修改调用方的输入切片
func TestVoxelFilterNoMutation(t *testing.T) {
	assert := assert.New(t)

	pos := []vec3d.T{{100, 200, 10}, {100.01, 200.01, 11}, {500, 700, 12}}
	orig := append([]vec3d.T{}, pos...)
	vg := newVoxelGrid(vec3d.T{10, 10, 10})
	_, err := vg.Filter(pos)
	assert.NoError(err)
	assert.Equal(orig, pos)
}

// 体素滤波必须用稀疏结构：稠密预分配在默认 filterSize 下是 21 GB 级分配。
// 这里用可控规模验证大间距网格不再分配 (xs+1)(ys+1)(zs+1) 规模的数组。
func TestVoxelFilterSparseMemory(t *testing.T) {
	assert := assert.New(t)

	// 两点距离极远，leafSize 极小：稠密实现会分配天文数字的体素
	pos := []vec3d.T{{0, 0, 0}, {1e6, 1e6, 1e6}}
	vg := newVoxelGrid(vec3d.T{1e-3, 1e-3, 1e-3})
	out, err := vg.Filter(pos)
	assert.NoError(err)
	assert.Len(out, 2)
}

func loadVerifyPoints(t *testing.T) []vec3d.T {
	t.Helper()
	f, err := os.Open("./test.json")
	assert.NoError(t, err)
	defer f.Close()
	b, err := io.ReadAll(f)
	assert.NoError(t, err)
	fcs, err := general.UnmarshalFeatureCollection(b)
	assert.NoError(t, err)
	it := NewKrigingInterpolator(Options{Input: fcs})
	pos, err := it.filter(it.extractPosion())
	assert.NoError(t, err)
	return pos
}
