### The general guidelines

- all tests are complete and pass the criteria
- the software meets all functional and non-functional requirements
- obtaining approval from all relevant stakeholders
- resolving all defects and issues identified during testing
- documentation is complete and up-to-date

### Do the following things before announcing the new release

- change "USING_DEV_VERSION" variable to False in __init__.py in the release branch.
- create a new folder for the documentation files of the new release in the "gh-pages" branch.
- make sure all the Jupyter notebooks work without errors.
- make sure all examples work without errors.
- make sure the DEV warnings are not showing, especially the first big one.
- check the version number is correct.
- make sure the docker images work correctly.
- check the online documentation for errors.
- run pytest against all testcases.
- check pip installation on Windows, macOS and Ubuntu.
- check conda installation on Windows, macOS and Ubuntu.
- check all gplately commands work as expected.
- check the list of plate models is documented correctly.
- check the pipy page https://pypi.org/project/gplately/, especially the images and links 