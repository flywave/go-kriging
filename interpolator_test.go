package kriging

import (
	"io/ioutil"
	"os"
	"testing"

	"github.com/flywave/go-cog"
	"github.com/flywave/go-geom/general"
)

func TestInterpolator(t *testing.T) {

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

	gtiff := cog.Read("./out.tif")

	if gtiff == nil {
		t.FailNow()
	}
}
