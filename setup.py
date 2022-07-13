from setuptools import find_packages, setup

setup(
	name = "WMAPLike",
	version = "0.1.0",
	description = "WMAP DR5 likelihood for cobaya",
	author = "Hidde Jense",
	author_email = "JenseH@cardiff.ac.uk",
	
	zip_safe = True,
	packages = find_packages(),
	python_requires = ">=3.7",
	install_requires = [
		"cobaya>=3.1.0",
		"astropy"
	],
	package_data = { "wmaplike" : [ "WMAPLike.yaml" ] }
)
