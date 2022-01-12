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

func (f *voxelGrid) Filter(pc []vec3d.T) ([]vec3d.T, error) {
	min, max, err := minMaxVec3(pc)
	if err != nil {
		return nil, err
	}

	size := max.Sub(&min)
	xs, ys, zs := int(size[0]/f.LeafSize[0]), int(size[1]/f.LeafSize[1]), int(size[2]/f.LeafSize[2])
	voxels := make([]voxel, (xs+1)*(ys+1)*(zs+1))

	var n int
	for i := range pc {
		p := pc[i].Sub(&min)
		x, y, z := int(p[0]/f.LeafSize[0]), int(p[1]/f.LeafSize[1]), int(p[2]/f.LeafSize[2])
		v := &voxels[x+xs*(y+ys*z)]
		if v.num == 0 {
			v.index = i
			n++
		}
		v.num++
		v.sum.Add(p)
		i++
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
