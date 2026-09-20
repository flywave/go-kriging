package kriging

import (
	"errors"
	"sort"

	vec3d "github.com/flywave/go3d/float64/vec3"
)

type voxelGrid struct {
	LeafSize vec3d.T
}

type voxelKey [3]int

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

// Filter 对点云做体素栅格降采样：每个体素输出一个点（多点体素取质心）。
// 体素用稀疏 map 存储而不是稠密数组——默认 filterSize 1024×1024×512 的稠密
// 预分配约 21.5 GB，三维数据必然 OOM。
// 输出按体素下标排序，结果确定；不会修改调用方的输入切片。
func (f *voxelGrid) Filter(pc []vec3d.T) ([]vec3d.T, error) {
	min, max, err := minMaxVec3(pc)
	if err != nil {
		return nil, err
	}

	size := max.Sub(&min)
	xs, ys := leafCount(size[0], f.LeafSize[0]), leafCount(size[1], f.LeafSize[1])
	zs := leafCount(size[2], f.LeafSize[2])

	voxels := make(map[voxelKey]*voxel)
	for idx := range pc {
		lx, ly, lz := pc[idx][0]-min[0], pc[idx][1]-min[1], pc[idx][2]-min[2]
		key := voxelKey{
			voxelIndex(lx, f.LeafSize[0], xs),
			voxelIndex(ly, f.LeafSize[1], ys),
			voxelIndex(lz, f.LeafSize[2], zs),
		}
		v, ok := voxels[key]
		if !ok {
			v = &voxel{index: idx}
			voxels[key] = v
		}
		v.num++
		v.sum[0] += lx
		v.sum[1] += ly
		v.sum[2] += lz
	}

	keys := make([]voxelKey, 0, len(voxels))
	for k := range voxels {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		a, b := keys[i], keys[j]
		if a[2] != b[2] {
			return a[2] < b[2]
		}
		if a[1] != b[1] {
			return a[1] < b[1]
		}
		return a[0] < b[0]
	})

	newPc := make([]vec3d.T, 0, len(keys))
	for _, k := range keys {
		v := voxels[k]
		if v.num > 1 {
			f := MulFloat(&v.sum, 1.0/float64(v.num))
			f.Add(&min)
			newPc = append(newPc, *f)
		} else {
			newPc = append(newPc, pc[v.index])
		}
	}

	return newPc, nil
}
