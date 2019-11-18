from distutils.core import setup

setup(
    name="scrap",
    version="0.3",
    packages=[
        "scrap",
        "scrap.tags",
        "scrap.chromium",
        "scrap.plotting",
        "scrap.dataset",
        "scrap.utils"
    ],
    entry_points={
        "console_scripts": [
            "scrap-preprocess=scrap.dataset.preprocess:preprocess",
            "scrap-init-ws=scrap.dataset.initialize:initialize"
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
