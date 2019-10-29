from distutils.core import setup

setup(
    name="scrap",
    version="0.2",
    packages=[
        "scrap",
        "scrap.tags",
        "scrap.chromium",
        "scrap.plotting"
    ],
    install_requires=[
        "pandas",
        "pepars",
        "sparsedat",
        "numpy",
        "plotly",
        "h5py"
    ]
)
