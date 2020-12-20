import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cartogrammetry",
    packages=setuptools.find_packages(),
    version="0.0.3",
    author="Ivo de Liefde",
    author_email="ivodeliefde@gmail.com",
    license="MIT",
    description="A Python package for creating block or circle style cartograms.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ivodeliefde/cartogrammetry",
    install_requires=[
        "numpy",
        "pandas",
        "pulp",
        "scikit-learn",
        "fiona",
        "pyproj" "shapely",
        "geopandas",
        "matplotlib",
        "tqdm",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
