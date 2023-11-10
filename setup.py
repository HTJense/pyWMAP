from setuptools import find_packages, setup

setup(
    name="WMAPLike",
    version="0.1.1",
    description="WMAP DR5 likelihood for cobaya",
    author="Hidde Jense",
    author_email="JenseH@cardiff.ac.uk",

    zip_safe=True,
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "cobaya>=3.3.0",
        "astropy"
    ],
    package_data={"wmaplike": ["*.yaml", "tests/*"]},
)
