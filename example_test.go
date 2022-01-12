package kriging

import (
	"fmt"
	"io/ioutil"
	"os"

	"github.com/flywave/go-geom/general"
	vec3d "github.com/flywave/go3d/float64/vec3"
)

var (
	pos = []vec3d.T{
		{117.99598607600996, 31.995986076009952, 45.986076009952846},
		{117.99622303211338, 31.99622303211338, 46.223032113384235},
		{118.00282145442502, 32.002821454425025, 52.821454425024626},
		{118.03919253247047, 32.03919253247046, 89.19253247046487},
		{117.98106280242764, 31.981062802427637, 31.062802427638776},
	}
)

func ExampleVariogram_Contour_Exponential() {
	ordinaryKriging := New(pos)
	ordinaryKriging.Train(Exponential, 0, 100)
	contourRectangle := ordinaryKriging.Contour(200, 200)
	fmt.Printf("%#v", contourRectangle.Contour[:10])
	// Output:
	// []float64{31.062802427639, 31.67443506838088, 32.27805611994289, 32.87380457354172, 33.461820447530236, 34.042244827188064, 34.61521990152406, 35.18088899653797, 35.73939660436969, 36.29088840795865}

}

func ExampleVariogram() {
	f, _ := os.Open("./test.json")

	json, _ := ioutil.ReadAll(f)

	fcs, _ := general.UnmarshalFeatureCollection(json)

	opts := Options{
		Input:  fcs,
		Output: "./out.tif",
	}

	ker := NewKrigingInterpolator(opts)

	npos := ker.extractPosion()

	min, max, _ := MinMaxVec3(npos)

	vg := NewVoxelGrid(vec3d.T{(max[0] - min[0]) / 500, (max[1] - min[1]) / 500, (max[2] - min[2]) / 30})

	npos1, _ := vg.Filter(npos)

	//npos1 = append(npos1, npos[len(npos1)-1])

	ordinaryKriging := New(npos1)
	ordinaryKriging.Train(Exponential, 0, 100)
	contourRectangle := ordinaryKriging.Contour(200, 200)
	fmt.Printf("%#v", contourRectangle.Contour[:10])
	// Output:
	// []float64{31.062802427639, 31.67443506838088, 32.27805611994289, 32.87380457354172, 33.461820447530236, 34.042244827188064, 34.61521990152406, 35.18088899653797, 35.73939660436969, 36.29088840795865}

}
