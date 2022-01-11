package kriging

import (
	"fmt"

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

func ExampleVariogram_Contour_Spherical() {
	ordinaryKriging := New(pos)
	ordinaryKriging.Train(Spherical, 0, 100)
	contourRectangle := ordinaryKriging.Contour(200, 200)
	fmt.Printf("%#v", contourRectangle.Contour[:10])
	// Output:
	// []float64{31.062802427637955, 31.35568613698695, 31.649070507173384, 31.942636980745615, 32.23631166432933, 32.53011270795399, 32.82405871070559, 33.11816872323164, 33.41246224930053, 33.70695924637588}
}

func ExampleVariogram_Contour_Gaussian() {
	ordinaryKriging := New(pos)
	ordinaryKriging.Train(Gaussian, 0, 100)
	contourRectangle := ordinaryKriging.Contour(200, 200)
	fmt.Printf("%#v", contourRectangle.Contour[:10])
	// Output:
	// []float64{31.062802438363697, 31.19418275387546, 31.32895545000455, 31.467084521380848, 31.60853364267539, 31.753266184588494, 31.901245229805934, 32.052433588854164, 32.206793815838, 32.36428822416533}
}
