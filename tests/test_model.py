from pytest import approx, raises

from datetime import datetime
from igrf_model import IGRFModel


def test_model_v13():
    model = IGRFModel.get(version=13)
    assert isinstance(model, IGRFModel)
    assert len(model._table) == (19 * 120 + 7 * 195 + 1)
    assert model._table[0] == -31543
    assert model._table[250] == -389


def test_caching():
    model1 = IGRFModel.get(version=13)
    model2 = IGRFModel.get(version=13)
    assert model1 is model2


def test_invalid_version():
    with raises(ValueError, match="unsupported model version"):
        IGRFModel.get(version=777777)


def test_evaluation():
    model = IGRFModel.get(version=13)

    # Budapest according to http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html,
    # on Dec 22, 2021
    vec = model.evaluate(47.498, 19.04, date=datetime(2021, 12, 22))
    assert vec.declination == approx(5.371, abs=1e-2)
    assert vec.inclination == approx(64.265, abs=1e-2)
    assert vec.total_intensity == approx(48978, abs=1)

    # Genova according to http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html,
    # 2000m above MSL, on Jun 30, 2002
    vec = model.evaluate(44.414, 8.942, alt=2000, date=datetime(2002, 6, 30))
    assert vec.declination == approx(0.658, abs=1e-2)
    assert vec.inclination == approx(60.457, abs=1e-2)
    assert vec.total_intensity == approx(46543, abs=1)

    # Tokyo according to http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html,
    # 1000m above MSL, on Jan 15, 1996
    vec = model.evaluate(35.652, 139.83, alt=1000, date=datetime(1996, 1, 15))
    assert vec.declination == approx(-6.844, abs=1e-2)
    assert vec.inclination == approx(48.902, abs=1e-2)
    assert vec.total_intensity == approx(46175, abs=1)

    # Buenos Aires according to http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html,
    # 1500m above MSL, on Jul 29, 1983
    vec = model.evaluate(-34.603, -58.381, alt=1500, date=datetime(1983, 7, 29))
    assert vec.declination == approx(-3.940, abs=1e-2)
    assert vec.inclination == approx(-34.068, abs=1e-2)
    assert vec.total_intensity == approx(24227, abs=1)


def test_evaluation_invalid_years():
    model = IGRFModel.get(version=13)

    with raises(ValueError, match="year"):
        model.evaluate(47.498, 19.04, date=datetime(1899, 12, 22))
    with raises(ValueError, match="year"):
        model.evaluate(47.498, 19.04, date=2050)


def test_evaluation_works_without_year():
    model = IGRFModel.get(version=13)
    now = datetime.now()
    result1 = model.evaluate(47.498, 19.04)
    result2 = model.evaluate(47.498, 19.04, date=now)
    assert result1 == approx(result2)
