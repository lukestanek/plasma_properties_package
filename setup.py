import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="plasma_properties_package", # Replace with your own username
    version="0.0.8",
    author="Lucas J. Stanek",
    author_email="staneklu@msu.edu",
    description="Compute Plasma Properties",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lukestanek/plasma_properties",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
