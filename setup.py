import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cartogrammetry",
    version="0.0.1",
    author="Ivo de Liefde",
    author_email="ivodeliefde@gmail.com",
    description="A Python package for creating block or circle style cartograms.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ivodeliefde/cartogrammetry",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
