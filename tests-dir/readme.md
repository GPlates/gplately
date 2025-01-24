There are four sub-folders in this "tests-dir". This folder is called "tests-dir" instead of "tests" to avoid potential conflicts with assorted things, such as Python packaging. 

- data

    This "data" folder contains the data files for the tests.

- debug_utils

    This "debug_utils" folder contains scripts which are used by developers during their debug process, such as reproducing bugs, performance tuning, experimenting concepts, etc. Some of the scripts should be refined and added to pytestcases folder when the time comes.

- pytestcases

    This "pytestcases" folder contains the formal test cases which run by Pytest automatically in GitHub Action.

- unittest

    This "unittest" folder contains scripts which need human interaction to verify the results. They are not a part of the formal automated test suites.