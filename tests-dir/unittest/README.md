This `unittest` folder has **nothing** to do with the [Python unittest model](https://docs.python.org/3/library/unittest.html). This folder contains scripts which need human interaction to verify the results. They are not a part of the formal automated test suites.

The Python scripts in this folder were designed to be executed manually by developers and the test results need to be verified by human, such as plotting maps and checking map details are correct. Other examples include test cases which take long time to run or download large files from Internet, such as making age grids. These scripts will not be executed automatically by default. Some of them should be refined and added to `pytestcases` folder in the future though. 

This folder also contains the scripts which are used by developers during their development process.

GPlately uses [Pytest](https://docs.pytest.org) to automate its test cases. See the `tests-dir/pytestcases` folder. Formal test suites should be placed in the `pytestcases` folder. 
