package kriging

import (
	"errors"

	vec3d "github.com/flywave/go3d/float64/vec3"
)

type voxelGrid struct {
	LeafSize vec3d.T
}

type voxel struct {
	sum   vec3d.T
	num   int
	index int
}

func newVoxelGrid(leafSize vec3d.T) *voxelGrid {
	vg := &voxelGrid{LeafSize: leafSize}
	return vg
}

func minMaxVec3(ra []vec3d.T) (vec3d.T, vec3d.T, error) {
	if len(ra) == 0 {
		return vec3d.T{}, vec3d.T{}, errors.New("no point")
	}
	min, max := ra[0], ra[0]
	for i := 1; i < len(ra); i++ {
		v := ra[i]
		for i := range v {
			if v[i] < min[i] {
				min[i] = v[i]
			}
			if v[i] > max[i] {
				max[i] = v[i]
			}
		}
	}
	return min, max, nil
}

func MulFloat(vec *vec3d.T, v float64) *vec3d.T {
	vec[0] *= v
	vec[1] *= v
	vec[2] *= v
	return vec
}

// leafCount 返回某个轴上的体素层数。
// 退化轴（该轴范围为 0，或 leafSize 非正）按单层处理：否则 size/leaf 会得到 NaN，
// 而 int(NaN) 的结果是实现相关的（amd64 上为 INT_MIN），会让下面的体素索引越界 panic。
// 扁平点集（所有点等高，即典型的“用 2D 点插值 Z”）必然走这条路。
func leafCount(size, leaf float64) int {
	if leaf <= 0 || size <= 0 {
		return 0
	}
	return int(size / leaf)
}

// voxelIndex 计算某个点在某轴上的体素下标，并夹在 [0, maxIndex] 内，
// 保证任何输入都不会算出越界下标
func voxelIndex(v, leaf float64, maxIndex int) int {
	if leaf <= 0 {
		return 0
	}
	i := int(v / leaf)
	if i < 0 {
		return 0
	}
	if i > maxIndex {
		return maxIndex
	}
	return i
}

func (f *voxelGrid) Filter(pc []vec3d.T) ([]vec3d.T, error) {
	min, max, err := minMaxVec3(pc)
	if err != nil {
		return nil, err
	}

	size := max.Sub(&min)
	xs, ys := leafCount(size[0], f.LeafSize[0]), leafCount(size[1], f.LeafSize[1])
	zs := leafCount(size[2], f.LeafSize[2])
	voxels := make([]voxel, (xs+1)*(ys+1)*(zs+1))

	var n int
	for idx := range pc {
		p := pc[idx].Sub(&min)
		x, y, z := voxelIndex(p[0], f.LeafSize[0], xs), voxelIndex(p[1], f.LeafSize[1], ys), voxelIndex(p[2], f.LeafSize[2], zs)
		v := &voxels[x+xs*(y+ys*z)]
		if v.num == 0 {
			v.index = idx
			n++
		}
		v.num++
		v.sum.Add(p)
	}

	newPc := make([]vec3d.T, 0, len(pc))
	for i := range voxels {
		v := &voxels[i]
		if n := v.num; n > 0 {
			if n > 1 {
				f := MulFloat(&v.sum, 1.0/float64(n))
				f.Add(&min)
				newPc = append(newPc, *f)
			} else {
				newPc = append(newPc, *pc[v.index].Add(&min))
			}
		}
	}

	return newPc, nil
}
