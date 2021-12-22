from pytest import approx

from datetime import datetime
from igrf_model import IGRFModel
from igrf_model.utils import parse_datetime_or_fractional_year


def test_date_property():
    now = datetime.now()
    model1 = IGRFModel.get().at(now)
    assert model1.date == parse_datetime_or_fractional_year(now)

    model2 = IGRFModel.get().now()
    assert model2.date == approx(model1.date, abs=1e-5)


def test_evaluation():
    model = IGRFModel.get(version=13)

    # Budapest according to http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html,
    # on Dec 22, 2021
    bound_model = model.at(datetime(2021, 12, 22))
    vec = bound_model.evaluate(47.498, 19.04)
    assert vec.declination == approx(5.371, abs=1e-2)
    assert vec.inclination == approx(64.265, abs=1e-2)
    assert vec.total_intensity == approx(48978, abs=1)

    # Genova according to http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html,
    # 2000m above MSL, on Jun 30, 2002
    bound_model = model.at(datetime(2002, 6, 30))
    vec = bound_model.evaluate(44.414, 8.942, alt=2000)
    assert vec.declination == approx(0.658, abs=1e-2)
    assert vec.inclination == approx(60.457, abs=1e-2)
    assert vec.total_intensity == approx(46543, abs=1)

    # Tokyo according to http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html,
    # 1000m above MSL, on Jan 15, 1996
    bound_model = model.at(datetime(1996, 1, 15))
    vec = bound_model.evaluate(35.652, 139.83, alt=1000)
    assert vec.declination == approx(-6.844, abs=1e-2)
    assert vec.inclination == approx(48.902, abs=1e-2)
    assert vec.total_intensity == approx(46175, abs=1)

    # Buenos Aires according to http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html,
    # 1500m above MSL, on Jul 29, 1983
    bound_model = model.at(datetime(1983, 7, 29))
    vec = bound_model.evaluate(-34.603, -58.381, alt=1500)
    assert vec.declination == approx(-3.940, abs=1e-2)
    assert vec.inclination == approx(-34.068, abs=1e-2)
    assert vec.total_intensity == approx(24227, abs=1)
