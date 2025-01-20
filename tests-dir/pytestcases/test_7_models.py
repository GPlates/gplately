import pytest, logging
from plate_model_manager import PlateModelManager
from gplately.commands.list_models import get_model_names

OFFICIALLY_SUPPORTED_MODEL_NAMES = [
    "Alfonso2024",
    "Muller2022",
    "Zahirovic2022",
    "Merdith2021",
    "Clennett2020",
    "Clennett2020_M2019",
    "Clennett2020_S2013",
    "Muller2019",
    "Young2018",
    "TorsvikCocks2017",
    "Matthews2016",
    "Matthews2016_pmag_ref",
    "Muller2016",
    "Scotese2016",
    "Zahirovic2016",
    "Gibbons2015",
    "Zahirovic2014",
    "Shephard2013",
    "Gurnis2012",
    "Seton2012",
    "Muller2008",
]

logger = logging.getLogger("TestLog")


@pytest.mark.parametrize("model_name", OFFICIALLY_SUPPORTED_MODEL_NAMES)
def test_model_availability(model_name):
    logger.info(model_name)
    m = PlateModelManager().get_model(model_name)
    assert m


def test_model_list():
    s1 = set([x.lower() for x in OFFICIALLY_SUPPORTED_MODEL_NAMES])
    s2 = set([y.lower() for y in get_model_names()])
    assert len(s1 - s2) == 0


def make_sure_zahirovic2022_and_muller2008_in_the_list():
    assert "Muller2008".lower() in [y.lower() for y in get_model_names()]
    assert "Zahirovic2022".lower() in [y.lower() for y in get_model_names()]
