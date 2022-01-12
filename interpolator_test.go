package kriging

import (
	"io/ioutil"
	"os"
	"testing"

	"github.com/flywave/go-cog"
	"github.com/flywave/go-geom/general"
)

func TestInterpolator(t *testing.T) {

	gtiff := cog.Read("./test.tif")

	if gtiff == nil {
		t.FailNow()
	}

	bbox := gtiff.GetBounds(0)

	if bbox.Min[0] == 0 {
		t.FailNow()
	}

	f, _ := os.Open("./test.json")

	json, _ := ioutil.ReadAll(f)

	fcs, _ := general.UnmarshalFeatureCollection(json)
	m := ModelType("spherical")

	opts := Options{
		Input:  fcs,
		Output: "./out.tif",
		Model:  &m,
	}

	ker := NewKrigingInterpolator(opts)

	err := ker.Process()

	if err != nil {
		t.FailNow()
	}

}
