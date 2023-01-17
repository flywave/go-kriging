package kriging

import (
	"io/ioutil"
	"os"
	"testing"

	"github.com/flywave/go-geom/general"
)

func TestInterpolator1(t *testing.T) {

	f, _ := os.Open("./test.json")

	json, _ := ioutil.ReadAll(f)

	fcs, _ := general.UnmarshalFeatureCollection(json)
	m := ModelType("spherical")

	opts := Options{
		Input:  fcs,
		Output: "./out1.tif",
		Model:  &m,
	}

	ker := NewKrigingInterpolator(opts)

	_, _, err := ker.Process()

	if err != nil {
		t.FailNow()
	}

}

func TestInterpolator(t *testing.T) {

	f, _ := os.Open("./02.geojson")

	json, _ := ioutil.ReadAll(f)

	fcs, _ := general.UnmarshalFeatureCollection(json)
	m := ModelType("spherical")

	bg := "./biguiyuan.tif"

	opts := Options{
		Input:      fcs,
		Output:     "./out21.tif",
		Model:      &m,
		Background: &bg,
	}

	ker := NewKrigingInterpolator(opts)

	_, _, err := ker.Process()

	if err != nil {
		t.FailNow()
	}

}
