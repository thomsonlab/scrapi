from setuptools import setup

setup(
    name="scrapi",
    version="0.5.1",
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
        "pandas>=0.25.0",
        "pepars",
        "sparsedat>=1.0.0alpha5",
        "numpy",
        "plotly",
        "h5py",
        "sklearn",
        "scipy",
        "statsmodels"
    ]
)
