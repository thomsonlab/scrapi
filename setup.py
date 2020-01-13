from setuptools import setup

setup(
    name="scrapi",
    version="0.4.0",
    packages=[
        "scrapi",
        "scrapi.tags",
        "scrapi.chromium",
        "scrapi.plotting",
        "scrapi.dataset",
        "scrapi.utils"
    ],
    entry_points={
        "console_scripts": [
            "scrap-preprocess=scrapi.dataset.preprocess:preprocess",
            "scrap-init-ws=scrapi.dataset.initialize:initialize"
        ]
    },
    install_requires=[
        "pandas",
        "pepars",
        "sparsedat",
        "numpy",
        "plotly",
        "h5py",
        "sklearn",
        "scipy",
        "statsmodels"
    ]
)
